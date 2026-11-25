#' Single EM step: E-step responsibilities + M-step update of A, mu, and all
#' non-parametric background/triggering components
#'
#'
#' @param events A `semihawkes_events` object (from `prepare_events()`,
#'   already passed through `add_adaptive_bandwidth()` / `add_bg_integral()`).
#' @param background A `semihawkes_spatial_base` object from the *previous*
#'   iteration's `background_spatial()` call (only its grid/marks are
#'   reused as the spatial domain; a fresh weighted surface is computed
#'   inside this step).
#' @param trend,weekly,monthly The *previous* iteration's
#'   `semihawkes_temporal_base` objects (weekly/monthly/trend).
#' @param temporal_trigger,spatial_trigger The *previous* iteration's
#'   triggering-kernel objects.
#' @param mu,A The *previous* iteration's scalar parameter estimates.
#' @param pairs A `trigger_pairs()` index (built once, reused every
#'   iteration -- it does not depend on the current parameter estimates).
#' @param cache_dir Directory populated by `trigger_spatial_marks()` (built
#'   once, reused every iteration -- depends only on `spatial_trigger`'s
#'   grid, which is fixed across iterations).
#' @param cell_size Passed to `background_spatial()`.
#' @param nprocs Passed to `trigger_evaluate()`: a whole number of parallel
#'   workers (`1`, the default, runs serially), `"auto"`, or a cluster from
#'   [parallel::makeCluster()]. A number or `"auto"` starts new workers at
#'   every call; to reuse them across steps, pass a cluster, as
#'   [fit_semihawkes()] does.
#' @param exclude_self Passed to `trigger_evaluate()`; default `TRUE`.
#' @param bgint,trig Optional. The results of `background_intensity()` and
#'   `trigger_evaluate()` for exactly the `background`, `trend`, `weekly`,
#'   `monthly` and triggering kernels passed in, as returned in `$bgint`
#'   and `$trig` by the previous call to `em_step()`. Supplying them saves
#'   evaluating the same quantities again at the start of this step; the
#'   result is the same. `NULL` (the default) evaluates them here. Values
#'   computed for any other inputs give wrong results.
#' @param marks Optional: the spatial marks in `cache_dir`, already read into
#'   a list with one logical matrix per event, as [fit_semihawkes()] does once
#'   per fit. Passed to [trigger_edge_correction()]; `NULL` (the default)
#'   reads them from `cache_dir` at this step.
#'
#' @return A list of class `"semihawkes_em_step"` with the updated
#'   `mu`, `A`, `llik`, `lambda_at_all`, `lambda_at_events`, `bgprobs`,
#'   the refreshed `background`, `trend`, `weekly`, `monthly`,
#'   `temporal_trigger`, `spatial_trigger` objects, and `bgint` and `trig`,
#'   the background intensity and triggering contributions evaluated with
#'   those refreshed components -- everything needed to feed straight into
#'   the next call to `em_step()`.
#' @examples
#' \donttest{
#' # A quick fit to one year of the bundled simulated data: a coarse grid, a
#' # coarse triggering-kernel grid and two EM iterations keep it to a few
#' # seconds. Real analyses need finer settings and more iterations.
#' cache_dir <- tempfile("spHawkes_marks_")
#' fit <- fit_semihawkes(
#'   sim_events, sim_region,
#'   date_col  = "date",
#'   end_date  = "2022-01-01",
#'   cell_size = 1000,
#'   cache_dir = cache_dir,
#'   trigger_spatial_args = list(x_range = 2000, x_step = 250),
#'   trigger_window_dist  = 2000,
#'   n_iter  = 2,
#'   verbose = FALSE,
#'   seed    = 1
#' )
#' # One more EM step from where the fit stopped
#' pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
#' step <- em_step(fit$events, fit$background, fit$trend, fit$weekly,
#'                 fit$monthly, fit$temporal_trigger, fit$spatial_trigger,
#'                 mu = fit$mu, A = fit$A, pairs = pairs,
#'                 cache_dir = cache_dir, cell_size = 1000)
#' step
#' }
#' @export
em_step <- function(events, background, trend, weekly, monthly,
                    temporal_trigger, spatial_trigger, mu, A,
                    pairs, cache_dir, cell_size, nprocs = 1, exclude_self = TRUE,
                    bgint = NULL, trig = NULL, marks = NULL) {

  ev <- events$events

  ## --- 1. E-step: current intensity decomposition & responsibilities -----
  # The previous step ends by evaluating both of these with the components
  # passed in here, so fit_semihawkes() hands them over rather than having
  # them evaluated a second time.
  if (is.null(bgint)) {
    bgint <- background_intensity(events, trend, weekly, monthly, background)
  } else if (length(bgint$bg_at_events) != nrow(ev)) {
    stop("`bgint` does not match `events`; leave it as NULL to evaluate it here.",
         call. = FALSE)
  }
  if (is.null(trig)) {
    trig <- trigger_evaluate(events, background, temporal_trigger, spatial_trigger,
                             exclude_self = exclude_self, nprocs = nprocs)
  } else if (length(trig$at_events) != nrow(ev)) {
    stop("`trig` does not match `events`; leave it as NULL to evaluate it here.",
         call. = FALSE)
  }

  lambda_at_events <- mu * bgint$bg_at_events + A * trig$at_events
  bgprobs <- mu * bgint$bg_at_events / lambda_at_events

  # Per-component background responsibility weights
  # used to re-estimate weekly/monthly/trend below.
  mub_events <- bgint$mub_events
  weW <- weekly$fun(bgint$marks_events$wday)   * mub_events / lambda_at_events
  weT <- trend$fun(ev$days)                    * mub_events / lambda_at_events
  weM <- monthly$fun(bgint$marks_events$month) * mub_events / lambda_at_events

  ## --- 2. M-step (non-parametric part): re-estimate temporal components --
  weekly_new  <- temporal_base_weekly(events,  resp_weight = weW)
  monthly_new <- temporal_base_monthly(events, resp_weight = weM)
  trend_new   <- temporal_base_trend(events,   resp_weight = weT)

  ## --- 3. M-step: re-estimate spatial background surface -----------------
  background_new <- background_spatial(events, cell_size = cell_size,
                                       resp_weight = bgprobs)

  ## --- 4. Refresh background intensity shape with the new components -----
  bgint_new <- background_intensity(events, trend_new, weekly_new, monthly_new,
                                    background_new)

  ## --- 5. M-step: re-estimate triggering kernels --------------------------
  edge_correction <- trigger_edge_correction(events, temporal_trigger, spatial_trigger,
                                             cache_dir, marks = marks)
  pair_weights <- trigger_pair_weights(events, pairs, A, lambda_at_events,
                                       temporal_trigger, spatial_trigger)

  temporal_trigger_new <- trigger_reestimate_temporal(
    events, pairs, pair_weights, edge_correction, temporal_trigger
  )
  spatial_trigger_new <- trigger_reestimate_spatial(
    events, pairs, pair_weights, edge_correction, spatial_trigger
  )

  ## --- 6. Refresh trigger evaluation with the new kernels -----------------
  trig_new <- trigger_evaluate(events, background_new, temporal_trigger_new,
                               spatial_trigger_new, exclude_self = exclude_self,
                               nprocs = nprocs)

  ## --- 7. M-step (parametric part): closed-form A, mu MLE update ---------
  n_events <- nrow(ev)
  A_new  <- (n_events - sum(bgprobs)) / trig_new$at_all
  mu_new <- (n_events - A_new * trig_new$at_all) / bgint_new$bg_at_all

  lambda_at_events_new <- mu_new * bgint_new$bg_at_events + A_new * trig_new$at_events
  lambda_at_all_new    <- mu_new * bgint_new$bg_at_all    + A_new * trig_new$at_all

  llik <- -(-sum(log(lambda_at_events_new)) + lambda_at_all_new)
  bgprobs_new <- mu_new * bgint_new$bg_at_events / lambda_at_events_new

  out <- list(
    mu = mu_new, A = A_new, llik = llik,
    lambda_at_all = lambda_at_all_new,
    lambda_at_events = lambda_at_events_new,
    bgprobs = bgprobs_new,
    background = background_new,
    trend = trend_new, weekly = weekly_new, monthly = monthly_new,
    temporal_trigger = temporal_trigger_new, spatial_trigger = spatial_trigger_new,
    bgint = bgint_new, trig = trig_new
  )
  class(out) <- "semihawkes_em_step"
  out
}

#' @export
print.semihawkes_em_step <- function(x, ...) {
  cat("<semihawkes_em_step>\n")
  cat(sprintf("  mu = %.6g, A = %.6g, llik = %.6g, lambda_at_all = %.6g\n",
              x$mu, x$A, x$llik, x$lambda_at_all))
  invisible(x)
}

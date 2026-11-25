#' Fit a semi-parametric self-exciting (Hawkes) point process by EM
#'
#' It runs Loop 0 (the
#' parametric/initial-guess step, previously a separate script) followed by
#' repeated calls to `em_step()` until convergence, and returns everything
#' in a single in-memory `"semihawkes_fit"` object -- with an *optional*
#' `checkpoint_dir` for anyone who still wants a per-iteration on-disk
#' record (e.g. for a long unattended run on a cluster).
#'
#' @param points An `sf` POINT layer of events (see `prepare_events()`).
#' @param region An `sf` POLYGON/MULTIPOLYGON study-area boundary.
#' @param date_col Name of the date column in `points`.
#' @param cell_size Spatial background grid cell size (required, no
#'   default -- see `background_spatial()`).
#' @param cache_dir Directory for `trigger_spatial_marks()`'s per-event
#'   cache (required, no default -- see `trigger_spatial_marks()`).
#' @param start_date,end_date,crs Passed to `prepare_events()`. The modelled
#'   study period always ends the day after the last event, whatever
#'   `end_date` is; see `prepare_events()`.
#' @param bandwidth_args A named list of extra arguments passed to
#'   `add_adaptive_bandwidth()` (e.g. `list(k = 5, min_dist = 150)`).
#' @param trigger_window_days,trigger_window_dist Passed to
#'   `trigger_pairs()`. `trigger_window_days` defaults to `90`.
#'   `trigger_window_dist` defaults to `NULL`, which uses the reach of the
#'   spatial kernel's grid: 1000 linear units (1 km for coordinates in
#'   metres) unless `trigger_spatial_args` sets another `x_range`. They cannot
#'   reach past the triggering kernels' grids: at most the temporal kernel's
#'   90-day cutoff and the spatial kernel's extent (see
#'   `trigger_spatial_args`). Choose the spatial reach with care: see
#'   *Choosing the triggering reach* below.
#' @param trigger_spatial_args A named list of extra arguments passed to
#'   [trigger_spatial_init()], which builds the grid the spatial triggering
#'   kernel is estimated on (e.g. `list(x_range = 1500, x_step = 100)`). The
#'   default grid reaches 1000 linear units from each event. Default
#'   `list()` (the function's own defaults).
#' @param n_iter Maximum number of EM iterations (the original `maxit`);
#'   default `50`.
#' @param tol Convergence tolerance on `abs(llik - llik_prev)` (the
#'   original `tolA`); default `0.001`.
#' @param tol_A Stop early if `A` drops below this (the original's
#'   `A < 0.001` stopping rule); default `0.001`.
#' @param seed Optional integer seed for the Monte Carlo edge-correction
#'   integral computed by [add_bg_integral()]. That integral is stochastic, so
#'   two runs with the same data differ slightly unless a seed is fixed here
#'   (or `set.seed()` is called beforehand). Default `NULL` (unseeded).
#' @param restore_rng Passed to [bg_integral()] when `seed` is supplied. If
#'   `TRUE` (the default), the global random-number state is put back
#'   afterwards, so the seed affects only the edge-correction integral. If
#'   `FALSE`, the generator is left in the state the seeded draw produces, as
#'   earlier versions of the package did.
#' @param nprocs Parallel workers for `trigger_spatial_marks()` and
#'   `trigger_evaluate()`, the slowest steps of a fit: a whole number (the
#'   default, `1`, runs serially), `"auto"` to choose from the machine, or a
#'   cluster from [parallel::makeCluster()], which is used and left running.
#'   Workers are started once for the whole fit and stopped at the end.
#'   `"auto"` uses one worker per physical core but one (two on a two-core
#'   machine), at most 2 while `R CMD check` runs, and runs serially for
#'   fewer than 1,000 events, where starting the workers costs more than they
#'   save. A number larger than the machine's cores is reduced to the number
#'   of cores, with a warning. The workers are socket clusters, which work on
#'   Windows, macOS and Linux, and the fit is the same, to the last bit,
#'   whatever their number.
#' @param checkpoint_dir Optional directory. If supplied, `em_step()`'s
#'   full output is `saveRDS()`-ed there after every iteration (the
#'   in-memory-object equivalent of the original's per-loop `.RData`
#'   dumps), so a crashed/killed long run can be resumed. Default `NULL`
#'   (no checkpointing).
#' @param verbose Print progress messages each iteration; default `TRUE`.
#' @param exclude_self Passed to `trigger_evaluate()`/`em_step()` for the
#'   whole EM run; default `TRUE`. Set to `FALSE` to reproduce the original
#'   script's `calculate_excite()` exactly.
#'
#' @return A list of class `"semihawkes_fit"` with elements:
#'   \describe{
#'     \item{events}{The `semihawkes_events` object used throughout.}
#'     \item{history}{A `data.frame` with one row per iteration:
#'       `iteration`, `mu`, `A`, `llik`, `lambda_at_all` -- the in-memory
#'       equivalent of the original script's `results_summary` built by
#'       re-reading every `Loop/LoopK/Results.RData` file.}
#'     \item{background, trend, weekly, monthly, temporal_trigger,
#'       spatial_trigger}{The final iteration's fitted components, each
#'       with a ready-to-use `$fun` evaluator.}
#'     \item{mu, A, llik, lambda_at_all, lambda_at_events, bgprobs}{The
#'       final iteration's scalar/vector estimates.}
#'     \item{n_iter}{Number of EM iterations actually run.}
#'     \item{stop_reason}{`"A_below_tol_A"`, `"llik_converged"`, or
#'       `"max_iter"`.}
#'   }
#'
#' @section Choosing the triggering reach:
#' `trigger_window_dist`, together with the `x_range` (and `y_range`) of the
#' spatial kernel's grid in `trigger_spatial_args`, sets how far from an event
#' the triggering kernel can reach. The default, 1000 linear units (1 km for
#' coordinates in metres), is a starting point rather than a recommendation
#' for every dataset, so consider the value carefully. Base it on how far one
#' event can plausibly trigger another in the process being studied, and
#' allow a few times that distance:
#' \itemize{
#'   \item A reach much wider than that lets the kernel's tail absorb
#'     clustering in the background, which overstates the branching ratio
#'     `A`. On simulated data with offspring a few hundred metres from their
#'     parent, a 4 km reach overstated `A` by about half.
#'   \item A reach that is too short cuts off genuine triggering, which is
#'     then attributed to the background and understates `A`.
#' }
#' Refitting with two or three reaches shows how sensitive `A` is to the
#' choice. Set the reach through the kernel's grid, for example
#' `trigger_spatial_args = list(x_range = 1500)`; `trigger_window_dist`
#' follows it unless it is set too.
#'
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
#' fit
#' summary(fit)
#' plot(fit, which = "history")
#' }
#' @export
fit_semihawkes <- function(points, region, date_col, cell_size, cache_dir,
                           start_date = NULL, end_date = NULL, crs = NULL,
                           bandwidth_args = list(),
                           trigger_window_days = 90, trigger_window_dist = NULL,
                           trigger_spatial_args = list(),
                           n_iter = 50, tol = 0.001, tol_A = 0.001, seed = NULL,
                           nprocs = 1, checkpoint_dir = NULL, verbose = TRUE,
                           exclude_self = TRUE, restore_rng = TRUE) {

  .check_flag(restore_rng, "restore_rng")
  .check_nprocs(nprocs)

  ## --- Part 1: data preparation (unchanged from the earlier port) --------
  events <- prepare_events(points, region, date_col,
                           start_date = start_date, end_date = end_date, crs = crs)
  events <- do.call(add_adaptive_bandwidth, c(list(x = events), bandwidth_args))
  events <- if (is.null(seed)) {
    add_bg_integral(events)
  } else {
    add_bg_integral(events, seed = seed, restore_rng = restore_rng)
  }

  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }

  ## --- Loop 0: initial parametric/non-parametric guesses ------------------
  if (isTRUE(verbose)) message("Loop 0: initial guesses")

  trend   <- temporal_base_trend(events)
  weekly  <- temporal_base_weekly(events)
  monthly <- temporal_base_monthly(events)
  background <- background_spatial(events, cell_size = cell_size)

  temporal_trigger <- trigger_temporal_init()
  spatial_trigger  <- do.call(trigger_spatial_init, trigger_spatial_args)

  # The temporal kernel is estimated by binning the time lags of the pairs
  # `trigger_pairs()` keeps onto the kernel's grid, which ends at its cutoff.
  # Longer lags cannot be binned: hist.weighted() used to report that from
  # inside `try()` and carry on, and the fit then failed with an unrelated
  # "length mismatch in convolution".
  temporal_cutoff <- max(temporal_trigger$base)
  if (!is.numeric(trigger_window_days) || length(trigger_window_days) != 1 ||
      is.na(trigger_window_days)) {
    stop("`trigger_window_days` must be a single number.", call. = FALSE)
  }
  if (trigger_window_days > temporal_cutoff) {
    stop(sprintf(
      paste0("`trigger_window_days` (%s) is longer than the temporal ",
             "triggering kernel's cutoff (%s days), so the time lags of ",
             "pairs beyond it cannot be binned. Use `trigger_window_days = %s` ",
             "or less."),
      format(trigger_window_days), format(temporal_cutoff), format(temporal_cutoff)),
      call. = FALSE)
  }

  # A kernel grid reaching further than the study region itself means most of
  # the grid holds offsets no pair of events can realise, so the kernel is
  # estimated from a small corner of it. Warn once here rather than per
  # iteration.
  kernel_reach <- min(max(spatial_trigger$x), max(spatial_trigger$y))
  if (is.null(trigger_window_dist)) {
    trigger_window_dist <- kernel_reach
  }

  # `trigger_pairs()` keeps pairs up to `trigger_window_dist` apart, and the
  # spatial kernel is estimated by binning those offsets onto the kernel grid.
  # If the window reaches past the grid, the out-of-range offsets cannot be
  # binned and the density matrix is built from a mismatched length -- wrong
  # numbers, with only a recycling warning to show for it.
  if (trigger_window_dist > kernel_reach) {
    stop(sprintf(
      paste0("`trigger_window_dist` (%.0f) is larger than the spatial ",
             "triggering kernel's extent (%.0f), so pair offsets beyond the ",
             "kernel grid cannot be binned. Set them to match, e.g. ",
             "`trigger_window_dist = %.0f` with ",
             "`trigger_spatial_args = list(x_range = %.0f)`."),
      trigger_window_dist, kernel_reach, kernel_reach, kernel_reach),
      call. = FALSE)
  }

  region_bbox <- sf::st_bbox(events$region)
  region_span <- max(region_bbox[["xmax"]] - region_bbox[["xmin"]],
                     region_bbox[["ymax"]] - region_bbox[["ymin"]])
  if (kernel_reach > region_span) {
    # A whole number of the current grid steps, as trigger_spatial_init()
    # requires of `x_range`.
    x_step <- spatial_trigger$x[2] - spatial_trigger$x[1]
    suggested <- x_step * max(1, round(signif(region_span / 2, 2) / x_step))
    warning(sprintf(
      paste0("The spatial triggering kernel reaches %.0f units from each ",
             "event, but the study region spans only %.0f, so most of that ",
             "grid holds offsets no pair of events can realise and the kernel ",
             "is estimated from a small corner of it. Consider ",
             "`trigger_spatial_args = list(x_range = %s)` together with ",
             "`trigger_window_dist = %s`."),
      kernel_reach, region_span, format(suggested), format(suggested)),
      call. = FALSE)
  }

  pairs <- trigger_pairs(events, window_days = trigger_window_days,
                         window_dist = trigger_window_dist)

  # One set of parallel workers for the whole fit: starting them takes about
  # as long as one trigger_evaluate() call on a few thousand events, which
  # would cancel out most of the gain if each call started its own.
  workers <- .open_workers(nprocs, nrow(events$events))
  if (workers$owned) on.exit(parallel::stopCluster(workers$cl), add = TRUE)
  if (!is.null(workers$cl)) {
    if (isTRUE(verbose)) message(sprintf("Using %d parallel workers", length(workers$cl)))
    nprocs <- workers$cl
  } else {
    nprocs <- 1L
  }

  trigger_spatial_marks(events, spatial_trigger, cache_dir = cache_dir, nprocs = nprocs)
  # Read once here rather than from disk at every iteration.
  marks <- .read_spatial_marks(cache_dir, nrow(events$events))

  bgint0 <- background_intensity(events, trend, weekly, monthly, background)
  trig0  <- trigger_evaluate(events, background, temporal_trigger, spatial_trigger,
                             exclude_self = exclude_self, nprocs = nprocs)

  n_events <- nrow(events$events)
  A  <- n_events / (2 * trig0$at_all)
  mu <- n_events / (2 * bgint0$bg_at_all)

  lambda_at_events <- mu * bgint0$bg_at_events + A * trig0$at_events
  lambda_at_all    <- mu * bgint0$bg_at_all    + A * trig0$at_all
  llik <- -(-sum(log(lambda_at_events)) + lambda_at_all)
  bgprobs <- mu * bgint0$bg_at_events / lambda_at_events

  history <- data.frame(iteration = 0L, mu = mu, A = A, llik = llik,
                        lambda_at_all = lambda_at_all)
  llik_prev <- llik
  stop_reason <- "max_iter"

  # Evaluated with the current components; each em_step() starts from these
  # and returns them re-evaluated with its updated components, so they are
  # computed once per iteration rather than twice.
  bgint <- bgint0
  trig  <- trig0

  n_iter_run <- 0L
  for (kk in seq_len(n_iter)) {
    n_iter_run <- kk
    if (isTRUE(verbose)) message(sprintf("Starting Loop: %d", kk))

    step <- em_step(events, background, trend, weekly, monthly,
                    temporal_trigger, spatial_trigger, mu, A,
                    pairs = pairs, cache_dir = cache_dir,
                    cell_size = cell_size, exclude_self = exclude_self,
                    nprocs = nprocs, bgint = bgint, trig = trig, marks = marks)

    mu <- step$mu; A <- step$A; llik <- step$llik
    lambda_at_all <- step$lambda_at_all
    lambda_at_events <- step$lambda_at_events
    bgprobs <- step$bgprobs
    background <- step$background
    trend <- step$trend; weekly <- step$weekly; monthly <- step$monthly
    temporal_trigger <- step$temporal_trigger
    spatial_trigger  <- step$spatial_trigger
    bgint <- step$bgint
    trig  <- step$trig

    history <- rbind(history, data.frame(
      iteration = kk, mu = mu, A = A, llik = llik, lambda_at_all = lambda_at_all
    ))

    if (isTRUE(verbose)) {
      message(sprintf("  mu = %.6g, A = %.6g, llik = %.6g", mu, A, llik))
    }

    if (!is.null(checkpoint_dir)) {
      saveRDS(step, file.path(checkpoint_dir, sprintf("step_%03d.rds", kk)))
    }

    if (A < tol_A) {
      stop_reason <- "A_below_tol_A"
      if (isTRUE(verbose)) message("A < tol_A: stopping.")
      break
    }
    if (!is.na(llik_prev) && abs(llik - llik_prev) < tol) {
      stop_reason <- "llik_converged"
      if (isTRUE(verbose)) {
        message(sprintf("|delta llik| = %.6g < %.6g: stopping.",
                        abs(llik - llik_prev), tol))
      }
      break
    }
    llik_prev <- llik
  }

  out <- list(
    events = events, history = history,
    background = background, trend = trend, weekly = weekly, monthly = monthly,
    temporal_trigger = temporal_trigger, spatial_trigger = spatial_trigger,
    mu = mu, A = A, llik = llik, lambda_at_all = lambda_at_all,
    lambda_at_events = lambda_at_events, bgprobs = bgprobs,
    n_iter = n_iter_run, stop_reason = stop_reason
  )
  class(out) <- "semihawkes_fit"
  out
}

#' @export
print.semihawkes_fit <- function(x, ...) {
  cat("<semihawkes_fit>\n")
  cat(sprintf("  %d iterations (stopped: %s)\n", x$n_iter, x$stop_reason))
  cat(sprintf("  mu = %.6g, A = %.6g, llik = %.6g\n", x$mu, x$A, x$llik))
  cat("  summary() for the background/triggering decomposition\n")
  invisible(x)
}

#' Summarize a fitted semi-parametric Hawkes model
#'
#' Reports what the `print()` method does not: the study window the model
#' was fitted over, the estimated branching structure (how the observed events
#' divide between background and triggered), the range of each fitted
#' background component, and how the EM run terminated.
#'
#' @section Interpreting `A`:
#' Because `g` and `h` are each normalized to integrate to 1, the branching
#' ratio `A` is the expected number of *direct* offspring generated by each
#' event. For a subcritical process (`A < 1`) it is also the long-run
#' proportion of events that are offspring rather than immigrants, so it can
#' be read directly as a triggered share. The empirical counterpart reported
#' alongside it is `1 - mean(bgprobs)`, the average fitted probability that an
#' observed event was triggered; the two agree when the model fits well and
#' the window is long relative to the triggering range.
#'
#' @param object A `semihawkes_fit` object from [fit_semihawkes()].
#' @param ... Currently unused.
#'
#' @return An object of class `"summary.semihawkes_fit"`, a list with
#'   elements `n_events`, `start_date`, `end_date`, `TT`, `crs`, `mu`, `A`,
#'   `llik`, `triggered_share`, `background_share`, `n_background`,
#'   `n_triggered`, `components` (a `data.frame` of per-component ranges),
#'   `n_iter`, `stop_reason` and `llik_change`. Printed by
#'   [print.summary.semihawkes_fit()].
#'
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
#' summary(fit)
#' }
#' @export
summary.semihawkes_fit <- function(object, ...) {
  bg <- object$bgprobs
  n <- length(bg)

  rng <- function(v) c(min(v, na.rm = TRUE), max(v, na.rm = TRUE))
  components <- data.frame(
    component = c("trend", "weekly", "monthly", "spatial background",
                  "trigger (temporal)", "trigger (spatial)"),
    min = c(rng(object$trend$basevalue)[1],
            rng(object$weekly$basevalue)[1],
            rng(object$monthly$basevalue)[1],
            rng(object$background$basevalue[object$background$marks > 0])[1],
            rng(object$temporal_trigger$basevalue)[1],
            rng(object$spatial_trigger$basevalue)[1]),
    max = c(rng(object$trend$basevalue)[2],
            rng(object$weekly$basevalue)[2],
            rng(object$monthly$basevalue)[2],
            rng(object$background$basevalue[object$background$marks > 0])[2],
            rng(object$temporal_trigger$basevalue)[2],
            rng(object$spatial_trigger$basevalue)[2]),
    stringsAsFactors = FALSE
  )

  h <- object$history
  llik_change <- if (nrow(h) > 1L) {
    h$llik[nrow(h)] - h$llik[nrow(h) - 1L]
  } else {
    NA_real_
  }

  out <- list(
    n_events         = n,
    start_date       = object$events$start_date,
    end_date         = object$events$end_date,
    TT               = object$events$TT,
    crs              = object$events$crs,
    mu               = object$mu,
    A                = object$A,
    llik             = object$llik,
    triggered_share  = 1 - mean(bg),
    background_share = mean(bg),
    n_background     = sum(bg),
    n_triggered      = n - sum(bg),
    components       = components,
    n_iter           = object$n_iter,
    stop_reason      = object$stop_reason,
    llik_change      = llik_change
  )
  class(out) <- "summary.semihawkes_fit"
  out
}

#' @param x A `summary.semihawkes_fit` object.
#' @rdname summary.semihawkes_fit
#' @export
print.summary.semihawkes_fit <- function(x, ...) {
  reason <- switch(x$stop_reason,
                   A_below_tol_A  = "A fell below tol_A",
                   llik_converged = "log-likelihood converged",
                   max_iter       = "reached n_iter without converging",
                   x$stop_reason)

  cat("Semi-parametric spatio-temporal Hawkes fit\n")
  cat("------------------------------------------------------------\n")
  cat(sprintf("Events        : %d\n", x$n_events))
  cat(sprintf("Study window  : %s to %s (%s%.4g days)\n",
              x$start_date, x$end_date,
              .modelled_to(x$start_date, x$end_date, x$TT), x$TT))
  cat(sprintf("CRS           : %s\n", .crs_label(x$crs)))

  cat("\nParameters\n")
  cat(sprintf("  mu (baseline rate)     : %.6g  per unit area per day\n", x$mu))
  cat(sprintf("  A  (branching ratio)   : %.4f  expected offspring per event\n", x$A))
  if (x$A >= 1) {
    cat("  NOTE: A >= 1 implies a supercritical (explosive) process.\n")
  }
  cat(sprintf("  log-likelihood         : %.6g\n", x$llik))

  cat("\nBranching structure (from fitted background probabilities)\n")
  cat(sprintf("  background : %8.1f events (%5.1f%%)\n",
              x$n_background, 100 * x$background_share))
  cat(sprintf("  triggered  : %8.1f events (%5.1f%%)\n",
              x$n_triggered, 100 * x$triggered_share))

  cat("\nFitted components (range over their support)\n")
  comp <- x$components
  for (i in seq_len(nrow(comp))) {
    cat(sprintf("  %-19s [%10.4g, %10.4g]\n",
                comp$component[i], comp$min[i], comp$max[i]))
  }

  cat("\nConvergence\n")
  cat(sprintf("  %d EM iteration%s, %s\n",
              x$n_iter, if (x$n_iter == 1L) "" else "s", reason))
  if (!is.na(x$llik_change)) {
    cat(sprintf("  final change in log-likelihood: %.6g\n", x$llik_change))
  }
  invisible(x)
}

#' Plot method for a fitted `semihawkes_fit` object
#'
#' Replaces the script's scattered `png()`/`plot()`/`ggplot()`/
#' `filled.contour()` calls and the final `patchwork` combination with a
#' single method covering weekly/monthly/trend, the spatial background, the
#' two triggering kernels, and the EM convergence history (log-likelihood,
#' mu, A, lambda).
#'
#' @param x A `semihawkes_fit` object.
#' @param which One of `"weekly"`, `"monthly"`, `"trend"`, `"background"`,
#'   `"trigger_temporal"`, `"trigger_spatial"`, `"history"`.
#' @param ... Currently unused.
#' @return `x`, invisibly. Called for its side effect of drawing the plot.
#' @importFrom graphics axis filled.contour par
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
#' plot(fit, which = "weekly")
#' plot(fit, which = "trigger_temporal")
#' plot(fit, which = "background")
#' }
#' @export
plot.semihawkes_fit <- function(x, which = "history", ...) {
  which <- match.arg(which, c("weekly", "monthly", "trend", "background",
                              "trigger_temporal", "trigger_spatial", "history"))
  switch(which,
         weekly = {
           plot(x$weekly$breaks, x$weekly$basevalue, type = "l",
                main = "Weekly background", xaxt = "n", xlab = "", ylab = "")
           axis(1, at = seq(0, 7, 1),
                labels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", ""))
         },
         monthly = {
           plot(x$monthly$breaks, x$monthly$basevalue, type = "l",
                main = "Monthly background", xaxt = "n", xlab = "", ylab = "")
           axis(1, at = 0:11, labels = month.abb)
         },
         trend = plot(x$trend$breaks, x$trend$basevalue, type = "l",
                      main = "Long-term trend", xlab = "Day", ylab = ""),

         background = {
           region <- x$events$region
           grid_cropped <- crop_background_spatial(x$background, region, crop = TRUE)
           grid_cropped <- .plot_safe_crs(grid_cropped)
           region <- .plot_safe_crs(region)

           p <- ggplot2::ggplot() +
             ggplot2::geom_sf(data = grid_cropped, ggplot2::aes(fill = .data$basevalue),
                              colour = NA) +
             ggplot2::geom_sf(data = region, fill = NA, colour = "black", linewidth = 0.7) +
             ggplot2::scale_fill_gradientn(
               colours = c("darkblue", "yellow", "orange", "red"),
               na.value = "white",
               name = expression(mu[b])
             ) +
             ggplot2::scale_x_continuous(breaks = NULL) +
             ggplot2::scale_y_continuous(breaks = NULL) +
             ggplot2::coord_sf(datum = NA) +
             ggplot2::labs(title = "Spatial background", x = "X", y = "Y") +
             ggplot2::theme_minimal() +
             ggplot2::theme(
               axis.text = ggplot2::element_blank(),
               axis.ticks = ggplot2::element_blank(),
               axis.title = ggplot2::element_blank(),
               legend.position = "right",
               plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
             )

           print(p)
         },

         trigger_temporal = plot(x$temporal_trigger$base, x$temporal_trigger$basevalue,
                                 type = "l", lwd = 2, main = "Temporal triggering kernel"),
         trigger_spatial = filled.contour(x$spatial_trigger$x, x$spatial_trigger$y,
                                          x$spatial_trigger$basevalue,
                                          main = "Spatial triggering kernel"),
         history = {
           op <- par(mfrow = c(2, 2))
           on.exit(par(op))
           plot(x$history$iteration, x$history$llik, type = "l", main = "Log-likelihood",
                xlab = "Iteration", ylab = "llik")
           plot(x$history$iteration, x$history$mu, type = "l", main = "mu",
                xlab = "Iteration", ylab = "mu")
           plot(x$history$iteration, x$history$A, type = "l", main = "A",
                xlab = "Iteration", ylab = "A")
           plot(x$history$iteration, x$history$lambda_at_all, type = "l",
                main = "Integral of lambda", xlab = "Iteration", ylab = "lambda")
         }
  )
  invisible(x)
}

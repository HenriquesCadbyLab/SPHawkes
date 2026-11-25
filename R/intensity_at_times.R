#' Total model intensity integrated over space
#'
#' Computes the fitted semi-parametric Hawkes intensity integrated over the
#' study region at a vector of times,
#' \deqn{\Lambda(t) = \int_W \lambda(s, t) \, ds.}
#'
#' The background contribution is
#' \deqn{\mu \mu_T(t) \mu_w(t) \mu_m(t)
#'       \int_W \mu_b(s) \, ds,}
#' while the triggering contribution is
#' \deqn{A \sum_{i:t_i < t} g(t - t_i)
#'       \int_W h(s - s_i) \, ds.}
#'
#' Spatial integrals use the same background grid and edge-correction marks
#' as [background_intensity()] and [trigger_evaluate()], so the results are
#' numerically consistent with `fit$lambda_at_all`.
#'
#' @param fit A fitted object of class `semihawkes_fit` returned by
#'   [fit_semihawkes()].
#' @param times Numeric vector of times, in days, at which to evaluate the
#'   spatially integrated intensity. Times must lie within the fitted
#'   observation window.
#' @param exclude_self Logical; if `TRUE` (default), an event does not
#'   contribute to its own triggering intensity.
#'
#' @return A numeric vector, of the same length as `times`, containing the
#'   model intensity integrated over the study region at each time.
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
#' # Expected events per day across the study region, on three days
#' intensity_at_times(fit, times = c(30, 180, 330))
#' }
#' @export
intensity_at_times <- function(fit,
                               times,
                               exclude_self = TRUE) {
  stopifnot(inherits(fit, "semihawkes_fit"))

  if (!is.numeric(times) || any(!is.finite(times))) {
    stop("`times` must be a finite numeric vector.", call. = FALSE)
  }

  if (!is.logical(exclude_self) || length(exclude_self) != 1 ||
      is.na(exclude_self)) {
    stop("`exclude_self` must be a single logical value.",
         call. = FALSE)
  }

  times <- sort(unique(as.numeric(times)))

  if (any(times < 0) || any(times > fit$events$TT)) {
    stop(
      "`times` must lie within the fitted observation window [0, ",
      fit$events$TT, "].",
      call. = FALSE
    )
  }

  events <- fit$events

  cell_area <- .grid_cell_area(fit$background, events$region)

  # Same spatial grid and in-region mask used by trigger_evaluate().
  grid_x <- fit$background$x %o% rep(1, length(fit$background$y))
  grid_y <- rep(1, length(fit$background$x)) %o% fit$background$y
  in_region <- fit$background$marks > 0

  # Precompute the spatial integral of each event's triggering kernel.
  # This does not depend on t.
  ev <- events$events
  spatial_integrals <- vapply(
    seq_len(nrow(ev)),
    function(i) {
      grid_spatial <- fit$spatial_trigger$fun(
        grid_x - ev$coorx[i],
        grid_y - ev$coory[i]
      )
      sum(in_region * grid_spatial) * cell_area
    },
    numeric(1)
  )

  # Calendar marks for the requested evaluation times.
  time_marks <- .calendar_marks(times, events$start_date)

  # Background intensity integrated over space, as a function of time.
  background_spatial_integral <-
    sum(fit$background$basevalue * fit$background$marks) * cell_area

  lambda_background <- fit$mu *
    fit$trend$fun(times) *
    fit$weekly$fun(time_marks$wday) *
    fit$monthly$fun(time_marks$month) *
    background_spatial_integral

  # Triggering intensity integrated over space, as a function of time.
  lambda_triggering <- vapply(
    times,
    function(t_k) {
      if (exclude_self) {
        parent_index <- which(ev$days < t_k)
      } else {
        parent_index <- which(ev$days <= t_k)
      }

      if (length(parent_index) == 0) {
        return(0)
      }

      time_lags <- t_k - ev$days[parent_index]
      temporal_weights <- fit$temporal_trigger$fun(time_lags)

      sum(temporal_weights * spatial_integrals[parent_index])
    },
    numeric(1)
  )

  fit$A * lambda_triggering + lambda_background
}

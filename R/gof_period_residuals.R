#' Spatial residuals by time period
#'
#' Computes spatial raw residuals separately for each calendar year or
#' calendar month in the fitted observation window.
#'
#' Each window runs from the first to the last day of its year or month,
#' cut to the fitted observation window at either end. Together the windows
#' cover the observation window exactly once, so the observed and the
#' expected counts both add up to their totals over the whole window.
#'
#' @param fit A fitted object of class `semihawkes_fit` returned by
#'   [fit_semihawkes()].
#' @param period Either `"year"` or `"month"`.
#'
#' @return A named list of `semihawkes_gof_spatial` objects, one for each
#'   calendar year (named like `"2021"`) or month (`"2021-03"`) that overlaps
#'   the fitted observation window, including any without events.
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
#' res <- gof_period_residuals(fit, period = "month")
#' plot(res)
#' }
#' @export
gof_period_residuals <- function(fit,
                                 period = c("year", "month")) {
  stopifnot(inherits(fit, "semihawkes_fit"))

  period <- match.arg(period)

  # Calendar periods, cut to the fitted observation window. Each window used
  # to run from the first to the last event date in its period instead, so
  # it missed the days before the first event and after the last one.
  first_day <- fit$events$start_date
  last_day <- fit$events$start_date + fit$events$TT - 1
  label_format <- if (period == "year") "%Y" else "%Y-%m"
  first_period <- as.Date(format(first_day, if (period == "year") "%Y-01-01" else "%Y-%m-01"))

  period_starts <- seq(first_period, last_day, by = period)
  period_ends <- seq(first_period, by = period,
                     length.out = length(period_starts) + 1)[-1] - 1

  window_starts <- pmax(period_starts, first_day)
  window_ends <- pmin(period_ends, last_day)

  out <- lapply(
    seq_along(period_starts),
    function(k) {
      gof_spatial_residuals(
        fit = fit,
        start_date = window_starts[k],
        end_date = window_ends[k]
      )
    }
  )

  names(out) <- format(period_starts, label_format)
  class(out) <- "semihawkes_gof_period"
  out
}

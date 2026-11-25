#' Time-transformed goodness-of-fit diagnostic
#'
#' Computes the model-implied cumulative expected count over time,
#' \deqn{\tau(t) = \int_0^t \Lambda(u) \, du,}
#' where \eqn{\Lambda(t)} is the fitted intensity integrated over the study
#' region, and compares it with the observed cumulative event count.
#'
#' Under a correctly specified temporal intensity, the transformed event
#' times are approximately a unit-rate Poisson process. The function also
#' returns pointwise KS-style acceptance bands for the observed cumulative
#' count.
#'
#' @param fit A fitted object of class `semihawkes_fit` returned by
#'   [fit_semihawkes()].
#' @param dt Temporal grid spacing, in days, used to approximate the
#'   time integral.
#' @param level Coverage level of the pointwise KS-style bands.
#'
#' @return An object of class `semihawkes_gof_time`: a list containing
#'   `grid`, `events`, `bands` and `level`.
#'
#' @importFrom stats approx qbeta stepfun
#' @importFrom utils tail
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
#' gof <- gof_time_transform(fit, dt = 1)
#' plot(gof)
#' }
#' @export
gof_time_transform <- function(fit,
                               dt = 1,
                               level = 0.95) {
  stopifnot(inherits(fit, "semihawkes_fit"))

  if (!is.numeric(dt) || length(dt) != 1 || !is.finite(dt) || dt <= 0) {
    stop("`dt` must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(level) || length(level) != 1 ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  events <- fit$events$events
  n_events <- nrow(events)
  total_time <- fit$events$TT

  time_grid <- seq(0, total_time, by = dt)

  if (tail(time_grid, 1) < total_time) {
    time_grid <- c(time_grid, total_time)
  }

  lambda_grid <- intensity_at_times(fit, times = time_grid)

  grid_intervals <- diff(c(0, time_grid))
  expected_cumulative <- cumsum(lambda_grid * grid_intervals)

  observed_cumulative <- stepfun(
    events$days,
    0:n_events
  )(time_grid)

  # Transform each observed event time by interpolating the cumulative
  # expected-count curve.
  transformed_event_times <- approx(
    x = time_grid,
    y = expected_cumulative,
    xout = events$days,
    method = "linear",
    rule = 2
  )$y

  counts <- 0:n_events
  alpha <- (1 - level) / 2

  bands <- data.frame(
    count = counts,
    lower = n_events * qbeta(alpha, counts + 1, n_events - counts + 1),
    upper = n_events * qbeta(1 - alpha, counts + 1, n_events - counts + 1)
  )

  out <- list(
    grid = data.frame(
      time = time_grid,
      transformed_time = expected_cumulative,
      observed_cumulative = observed_cumulative,
      expected_cumulative = expected_cumulative
    ),
    events = data.frame(
      days = events$days,
      transformed_time = transformed_event_times
    ),
    bands = bands,
    level = level,
    dt = dt
  )

  class(out) <- "semihawkes_gof_time"
  out
}

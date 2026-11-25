#' Integrate gridded intensity over one time period
#'
#' Numerically integrates background and triggering intensity snapshots
#' returned by background_trigger_grid(), using the trapezoidal rule.
#'
#' Integration is performed independently at each spatial grid location.
#' Actual time differences are used, so equally spaced times are not
#' required.
#'
#' Both period boundaries must be available in grid$days. No temporal
#' interpolation or extrapolation is performed.
#'
#' This function integrates sampled intensity surfaces. It does not
#' reevaluate the fitted model or implement event-aware quadrature.
#' Accuracy should be assessed by refining the evaluation times,
#' particularly near triggering discontinuities.
#'
#' @param grid A list returned by background_trigger_grid().
#' @param start Initial model time, in days.
#' @param end Final model time, in days. Must be greater than start.
#'
#' @return A list with days, x, y, inside and snapshots, compatible with
#'   plot_intensity_snapshot(). The single snapshot contains temporally
#'   integrated bg and tr matrices. Its time identifier is start.
#'   Additional fields statistic, period and integration describe
#'   the aggregation.
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
#' grid <- background_trigger_grid(fit, sim_region, cell_size = 500,
#'                                 days = 0:30)
#' january <- integrate_intensity_period(grid, start = 0, end = 30)
#' str(january, max.level = 1)
#' }
#' @export
integrate_intensity_period <- function(grid, start, end) {

  required <- c("days", "x", "y", "inside", "snapshots")

  if (!is.list(grid) || !all(required %in% names(grid))) {
    stop(
      "'grid' must be an output of background_trigger_grid().",
      call. = FALSE
    )
  }

  if (!is.null(grid$statistic) &&
      !identical(grid$statistic, "snapshot")) {
    stop("Only instantaneous snapshots can be integrated.",
         call. = FALSE)
  }

  valid_time <- function(z) {
    is.numeric(z) &&
      length(z) == 1L &&
      is.finite(z)
  }

  if (!valid_time(start) || !valid_time(end) || end <= start) {
    stop(
      "'start' and 'end' must be finite times with 'end' > 'start'.",
      call. = FALSE
    )
  }

  if (!is.numeric(grid$days) || length(grid$days) < 2L ||
      any(!is.finite(grid$days)) || anyDuplicated(grid$days) ||
      length(grid$snapshots) != length(grid$days)) {
    stop("The grid must contain at least two distinct valid times.",
         call. = FALSE)
  }

  valid_axis <- function(z) {
    is.numeric(z) &&
      length(z) > 0L &&
      all(is.finite(z)) &&
      (length(z) == 1L || all(diff(z) > 0))
  }

  if (!valid_axis(grid$x) || !valid_axis(grid$y)) {
    stop("Grid coordinates must be finite and strictly increasing.",
         call. = FALSE)
  }

  expected_dim <- c(length(grid$x), length(grid$y))
  n_locations <- prod(expected_dim)

  if (!is.logical(grid$inside) ||
      length(grid$inside) != n_locations ||
      anyNA(grid$inside) || !any(grid$inside)) {
    stop("The grid must contain a valid non-empty inside-region mask.",
         call. = FALSE)
  }

  if (!all(c(start, end) %in% grid$days)) {
    stop(
      "Both period boundaries must be available in 'grid$days'. ",
      "Include them when evaluating background_trigger_grid().",
      call. = FALSE
    )
  }

  selected <- which(
    grid$days >= start & grid$days <= end
  )

  selected <- selected[order(grid$days[selected])]
  times <- grid$days[selected]
  dt <- diff(times)

  extract_component <- function(i, component) {
    z <- grid$snapshots[[i]][[component]]

    if (!is.matrix(z) || !is.numeric(z) ||
        !identical(dim(z), expected_dim)) {
      stop("A snapshot contains invalid intensity matrices.",
           call. = FALSE)
    }

    inside_values <- as.vector(z)[grid$inside]

    if (any(!is.finite(inside_values)) ||
        any(inside_values < 0)) {
      stop(
        "Inside-region intensity values must be finite and non-negative.",
        call. = FALSE
      )
    }

    # Outside-region locations do not contribute to the output.
    z[!matrix(grid$inside, nrow = expected_dim[1L],
              ncol = expected_dim[2L])] <- 0

    z
  }

  integrate_component <- function(component) {
    result <- matrix(
      0,
      nrow = expected_dim[1L],
      ncol = expected_dim[2L]
    )

    previous <- extract_component(selected[1L], component)

    for (k in seq_along(dt)) {
      current <- extract_component(selected[k + 1L], component)

      result <- result +
        dt[k] * (previous / 2 + current / 2)

      previous <- current
    }

    if (any(!is.finite(result))) {
      stop("Temporal integration produced non-finite values.",
           call. = FALSE)
    }

    result
  }

  integrated_background <- integrate_component("bg")
  integrated_triggering <- integrate_component("tr")

  list(
    days = start,
    x = grid$x,
    y = grid$y,
    inside = grid$inside,
    snapshots = list(
      list(
        bg = integrated_background,
        tr = integrated_triggering
      )
    ),
    statistic = "time_integral",
    period = c(start = start, end = end),
    integration = list(
      method = "trapezoid",
      times = times,
      time_unit = "days"
    )
  )
}

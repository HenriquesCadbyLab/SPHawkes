#' Spatial raw residuals for a fitted model
#'
#' Computes observed minus expected event counts on the spatial background
#' grid used by the fitted semi-parametric Hawkes model.
#'
#' For each grid cell \eqn{B}, the raw residual is
#' \deqn{R(B) = N(B) - \mathbb{E}[N(B)],}
#' where \eqn{N(B)} is the observed number of events in the cell and
#' \eqn{\mathbb{E}[N(B)]} is the model-implied expected count.
#'
#' The function also returns the relative residual,
#' \deqn{R_{\mathrm{rel}}(B) = \frac{N(B) - \mathbb{E}[N(B)]}
#'   {\mathbb{E}[N(B)]},}
#' and the Poisson-standardized residual,
#' \deqn{R_{\mathrm{std}}(B) = \frac{N(B) - \mathbb{E}[N(B)]}
#'   {\sqrt{\mathbb{E}[N(B)]}}.}
#'
#' @param fit A fitted object of class `semihawkes_fit` returned by
#'   [fit_semihawkes()].
#' @param start_date Optional first day of the evaluation window. If `NULL`,
#'   the start of the fitted observation window is used.
#' @param end_date Optional last day of the evaluation window, included in
#'   it: events on that day are counted, and the expected count runs to the
#'   end of the day. If `NULL`, the last day of the fitted observation window
#'   is used.
#'
#' @return An object of class `semihawkes_gof_spatial`: a list containing
#'   `cells`, `region`, `start_date` and `end_date` (the first and last days
#'   of the window), `start_day` and `end_day` (the window as the interval
#'   `[start_day, end_day)` in days since the start of the study) and
#'   `window_length` (`end_day - start_day`), and `cell_size`, the side of
#'   the background grid's cells. The `cells` data frame contains `x`, `y`,
#'   `expected`, `observed`, `residual`, `relative_residual`,
#'   `standardized_residual`.
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
#' res <- gof_spatial_residuals(fit)
#' plot(res)
#' }
#' @export
gof_spatial_residuals <- function(fit,
                                  start_date = NULL,
                                  end_date = NULL) {
  stopifnot(inherits(fit, "semihawkes_fit"))

  events <- fit$events
  ev <- events$events

  # Event `days` are whole days since the start of the study, and day `d`
  # covers the interval [d, d + 1). `start_date` and `end_date` are the first
  # and last days of the window, both included, so the window is the interval
  # [start_day, end_day) with `end_day` one past the last day. Events on
  # `end_date` used to be counted as observed while the expected count stopped
  # at the start of that day.
  start_day <- 0
  last_day <- events$TT - 1

  if (!is.null(start_date)) {
    start_day <- as.numeric(.local_date(start_date) - events$start_date)
  }

  if (!is.null(end_date)) {
    last_day <- as.numeric(.local_date(end_date) - events$start_date)
  }

  if (start_day < 0 || last_day > events$TT - 1 || start_day > last_day) {
    stop(
      "`start_date` and `end_date` must define a non-empty window within ",
      "the fitted observation window.",
      call. = FALSE
    )
  }

  end_day <- last_day + 1
  window_length <- end_day - start_day

  cell_area <- .grid_cell_area(fit$background, events$region)

  time_marks <- events$time.marks
  time_marks_window <- time_marks[
    time_marks >= start_day & time_marks <= end_day
  ]

  calendar_marks <- .calendar_marks(time_marks_window, events$start_date)

  # Background expected count in each grid cell.
  background_temporal_integral <-
    mean(
      fit$trend$fun(time_marks_window) *
        fit$weekly$fun(calendar_marks$wday) *
        fit$monthly$fun(calendar_marks$month)
    ) * window_length

  expected_background <- fit$mu *
    background_temporal_integral *
    fit$background$basevalue *
    fit$background$marks *
    cell_area

  # Triggering expected count in each grid cell.
  grid_x <- fit$background$x %o% rep(1, length(fit$background$y))
  grid_y <- rep(1, length(fit$background$x)) %o% fit$background$y
  in_region <- fit$background$marks > 0

  # Only events before the end of the window can trigger events in it, and
  # only those within the temporal kernel's reach of its start: the kernel is
  # zero beyond its cutoff, so earlier events add exactly zero. The day of
  # margin keeps rounding in the time marks from dropping an event right at
  # the cutoff.
  temporal_reach <- max(fit$temporal_trigger$base)
  parent_events <- ev[ev$days < end_day &
                        ev$days >= start_day - temporal_reach - 1, ]

  expected_triggering <- matrix(
    0,
    nrow = nrow(grid_x),
    ncol = ncol(grid_y)
  )

  # Each event's triggering is evaluated only on the cells within the spatial
  # kernel's grid around it: the kernel is zero outside that grid, so the
  # other cells would add exactly zero.
  kernel_x <- range(fit$spatial_trigger$x)
  kernel_y <- range(fit$spatial_trigger$y)

  if (nrow(parent_events) > 0) {
    for (i in seq_len(nrow(parent_events))) {
      temporal_integral <- mean(
        fit$temporal_trigger$fun(time_marks_window - parent_events$days[i])
      ) * window_length
      if (temporal_integral == 0) next

      dx <- fit$background$x - parent_events$coorx[i]
      dy <- fit$background$y - parent_events$coory[i]
      ix <- which(dx >= kernel_x[1] & dx <= kernel_x[2])
      iy <- which(dy >= kernel_y[1] & dy <= kernel_y[2])
      if (length(ix) == 0 || length(iy) == 0) next

      spatial_values <- fit$spatial_trigger$fun(
        grid_x[ix, iy, drop = FALSE] - parent_events$coorx[i],
        grid_y[ix, iy, drop = FALSE] - parent_events$coory[i]
      )

      expected_triggering[ix, iy] <- expected_triggering[ix, iy] +
        temporal_integral * in_region[ix, iy, drop = FALSE] * spatial_values * cell_area
    }
  }

  expected <- expected_background + fit$A * expected_triggering

  # Map each event to its nearest background-grid cell.
  x_index <- round(
    (ev$coorx - min(fit$background$x)) / fit$background$cell_size
  ) + 1

  y_index <- round(
    (ev$coory - min(fit$background$y)) / fit$background$cell_size
  ) + 1

  valid_index <- x_index >= 1 & x_index <= nrow(grid_x) &
    y_index >= 1 & y_index <= ncol(grid_y) &
    ev$days >= start_day & ev$days < end_day

  observed <- matrix(
    0,
    nrow = nrow(grid_x),
    ncol = ncol(grid_y)
  )

  if (any(valid_index)) {
    for (i in which(valid_index)) {
      observed[x_index[i], y_index[i]] <-
        observed[x_index[i], y_index[i]] + 1
    }
  }

  residual <- observed - expected
  relative_residual <- residual / expected
  standardized_residual <- residual / sqrt(expected)

  cells <- data.frame(
    x = rep(fit$background$x, times = length(fit$background$y)),
    y = rep(fit$background$y, each = length(fit$background$x)),
    expected = as.vector(expected),
    observed = as.vector(observed),
    residual = as.vector(residual),
    relative_residual = as.vector(relative_residual),
    standardized_residual = as.vector(standardized_residual),
    in_region = as.vector(in_region)
  )

  cells <- cells[cells$in_region, ]
  cells$in_region <- NULL

  out <- list(
    cells = cells,
    region = events$region,
    start_date = events$start_date + start_day,
    end_date = events$start_date + last_day,
    start_day = start_day,
    end_day = end_day,
    window_length = window_length,
    cell_size = fit$background$cell_size
  )

  class(out) <- "semihawkes_gof_spatial"
  out
}

#' Plot a time-transformed goodness-of-fit diagnostic
#'
#' Plots the observed cumulative event count against the model-implied
#' transformed time, together with pointwise KS-style acceptance bands.
#'
#' @param x An object of class `semihawkes_gof_time` returned by
#'   [gof_time_transform()].
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A `ggplot` object.
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
#' plot(gof_time_transform(fit))
#' }
#' @export
plot.semihawkes_gof_time <- function(x, ...) {
  stopifnot(inherits(x, "semihawkes_gof_time"))

  n_events <- max(x$bands$count)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = x$bands,
      mapping = ggplot2::aes(
        x = .data$count,
        ymin = .data$lower,
        ymax = .data$upper
      ),
      fill = "grey80",
      alpha = 0.5
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    ggplot2::geom_step(
      data = x$grid,
      mapping = ggplot2::aes(
        x = .data$transformed_time,
        y = .data$observed_cumulative
      ),
      colour = "#440154",
      linewidth = 0.6
    ) +
    ggplot2::coord_equal(
      xlim = c(0, n_events),
      ylim = c(0, n_events),
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "Transformed time",
      y = "Cumulative event count",
      title = "Time-rescaling goodness-of-fit diagnostic",
      subtitle = paste0(
        "Pointwise KS-style bands at level ",
        sprintf("%.0f%%", 100 * x$level)
      )
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

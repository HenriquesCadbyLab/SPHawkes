#' Plot period-wise spatial residuals
#'
#' Plots spatial residuals for each time period using a common colour scale,
#' so residual maps are directly comparable.
#'
#' Grid cells are calculated on the fitted background grid. For display,
#' cells intersecting the study boundary are clipped to the region polygon;
#' this affects visualisation only, not residual calculation.
#'
#' @param x An object of class `semihawkes_gof_period` returned by
#'   [gof_period_residuals()].
#' @param type Type of residual to plot: `"raw"` for observed minus expected,
#'   `"relative"` for `(observed - expected) / expected`, or
#'   `"standardized"` for `(observed - expected) / sqrt(expected)`.
#' @param limits Optional numeric vector of length two giving common
#'   colour-scale limits. If `NULL`, limits are chosen from the 2nd and 98th
#'   percentiles of all selected residuals.
#' @param ncol Number of columns in the faceted plot.
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A `ggplot` object.
#'
#' @importFrom stats quantile
#' @importFrom rlang .data
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
#' plot(gof_period_residuals(fit, period = "month"))
#' }
#' @export
plot.semihawkes_gof_period <- function(x,
                                       type = c(
                                         "raw",
                                         "relative",
                                         "standardized"
                                       ),
                                       limits = NULL,
                                       ncol = 3,
                                       ...) {
  stopifnot(inherits(x, "semihawkes_gof_period"))

  if (length(x) == 0) {
    stop("`x` contains no residual periods.", call. = FALSE)
  }

  type <- match.arg(type)

  residual_column <- switch(
    type,
    raw = "residual",
    relative = "relative_residual",
    standardized = "standardized_residual"
  )

  residual_data <- lapply(
    names(x),
    function(period) {
      cells <- x[[period]]$cells

      data.frame(
        x = cells$x,
        y = cells$y,
        residual = cells[[residual_column]],
        period = period
      )
    }
  )

  residual_data <- do.call(rbind, residual_data)

  residual_data$period <- factor(
    residual_data$period,
    levels = names(x)
  )

  cell_dims <- .residual_cell_size(x[[1]]$cell_size, residual_data$x, residual_data$y)
  cell_width <- cell_dims[1]
  cell_height <- cell_dims[2]

  cell_polygons <- lapply(seq_len(nrow(residual_data)), function(i) {
    xi <- residual_data$x[i]
    yi <- residual_data$y[i]

    sf::st_polygon(list(matrix(
      c(
        xi - cell_width / 2, yi - cell_height / 2,
        xi + cell_width / 2, yi - cell_height / 2,
        xi + cell_width / 2, yi + cell_height / 2,
        xi - cell_width / 2, yi + cell_height / 2,
        xi - cell_width / 2, yi - cell_height / 2
      ),
      ncol = 2,
      byrow = TRUE
    )))
  })

  region <- x[[1]]$region

  residual_sf <- sf::st_sf(
    residual_data,
    geometry = sf::st_sfc(
      cell_polygons,
      crs = sf::st_crs(region)
    )
  )

  residual_sf <- suppressWarnings(sf::st_intersection(residual_sf, region))
  residual_sf <- .plot_safe_crs(residual_sf)
  region <- .plot_safe_crs(region)

  if (is.null(limits)) {
    limits <- stats::quantile(
      residual_sf$residual,
      probs = c(0.02, 0.98),
      na.rm = TRUE
    )
  }

  residual_label <- switch(
    type,
    raw = "Residual",
    relative = "Relative residual",
    standardized = "Standardized residual"
  )

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = residual_sf,
      mapping = ggplot2::aes(fill = .data$residual),
      colour = NA
    ) +
    ggplot2::geom_sf(
      data = region,
      inherit.aes = FALSE,
      fill = NA,
      colour = "black",
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#1b7837",
      midpoint = 0,
      limits = limits,
      oob = scales::squish,
      na.value = NA,
      name = residual_label,
      guide = ggplot2::guide_colourbar(
        direction = "horizontal",
        barwidth = grid::unit(7, "cm"),
        barheight = grid::unit(0.35, "cm"),
        title.position = "top"
      )
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::facet_wrap(
      ~ period,
      ncol = ncol
    ) +
    ggplot2::labs(
      title = "Spatial residuals by period",
      subtitle = paste0("Residual type: ", type),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(hjust = 0.5),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

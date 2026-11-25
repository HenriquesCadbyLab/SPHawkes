#' Plot spatial residuals
#'
#' Plots observed minus expected event counts, or a relative or
#' Poisson-standardized version of that residual, on the spatial background
#' grid used by a fitted semi-parametric Hawkes model.
#'
#' Grid cells are calculated on the fitted background grid. For display,
#' cells intersecting the study boundary are clipped to the region polygon;
#' this affects visualisation only, not residual calculation.
#'
#' @param x An object of class `semihawkes_gof_spatial` returned by
#'   [gof_spatial_residuals()].
#' @param type Type of residual to plot: `"raw"` for observed minus expected,
#'   `"relative"` for `(observed - expected) / expected`, or
#'   `"standardized"` for `(observed - expected) / sqrt(expected)`.
#' @param limits Optional numeric vector of length two giving the colour-scale
#'   limits. If `NULL`, limits are chosen from the 2nd and 98th residual
#'   percentiles.
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
#' res <- gof_spatial_residuals(fit)
#' plot(res)
#' plot(res, type = "standardized")
#' }
#' @export
plot.semihawkes_gof_spatial <- function(x,
                                        type = c(
                                          "raw",
                                          "relative",
                                          "standardized"
                                        ),
                                        limits = NULL,
                                        ...) {
  stopifnot(inherits(x, "semihawkes_gof_spatial"))

  type <- match.arg(type)

  residual_column <- switch(
    type,
    raw = "residual",
    relative = "relative_residual",
    standardized = "standardized_residual"
  )

  plot_data <- x$cells
  plot_data$residual <- plot_data[[residual_column]]

  cell_dims <- .residual_cell_size(x$cell_size, plot_data$x, plot_data$y)
  cell_width <- cell_dims[1]
  cell_height <- cell_dims[2]

  cell_polygons <- lapply(seq_len(nrow(plot_data)), function(i) {
    xi <- plot_data$x[i]
    yi <- plot_data$y[i]

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

  cell_sf <- sf::st_sf(
    plot_data,
    geometry = sf::st_sfc(
      cell_polygons,
      crs = sf::st_crs(x$region)
    )
  )

  cell_sf <- suppressWarnings(sf::st_intersection(cell_sf, x$region))
  cell_sf <- .plot_safe_crs(cell_sf)
  plot_region <- .plot_safe_crs(x$region)

  if (is.null(limits)) {
    limits <- stats::quantile(
      cell_sf$residual,
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
      data = cell_sf,
      mapping = ggplot2::aes(fill = .data$residual),
      colour = NA
    ) +
    ggplot2::geom_sf(
      data = plot_region,
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
    ggplot2::labs(
      title = "Spatial residuals",
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
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# Width and height of the residual cells for the plot methods: the background
# grid's `cell_size`, which gof_spatial_residuals() stores. The plots used to
# take the gap between the first two x (and y) values of the cells in the
# region instead, which is two or more cells wide when the region has no
# cell in the second column (or row) of the grid. Residual objects made
# before `cell_size` was stored fall back to the smallest gap.
.residual_cell_size <- function(cell_size, x, y) {
  if (!is.null(cell_size)) {
    return(c(cell_size, cell_size))
  }
  smallest_gap <- function(v) {
    gaps <- diff(sort(unique(v)))
    if (length(gaps) > 0) min(gaps) else NA_real_
  }
  c(smallest_gap(x), smallest_gap(y))
}

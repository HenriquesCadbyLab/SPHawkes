#' Plot one gridded intensity snapshot
#'
#' Plots a previously evaluated background, triggering, or total intensity
#' snapshot returned by background_trigger_grid().
#'
#' Matrix rows correspond to x coordinates and matrix columns correspond
#' to y coordinates. Matrices are flattened in the same order as
#' expand.grid(x = grid$x, y = grid$y); no transpose is required.
#'
#' Only grid locations marked as inside the study region are displayed.
#' Tiles are centred on evaluation coordinates and are not clipped to the
#' study-region boundary. Values are pointwise intensities, not spatial
#' integrals over cells.
#'
#' @param grid A list returned by background_trigger_grid().
#' @param time A single evaluated model time, in days. This selects a
#'   value in grid$days, not a snapshot position.
#' @param component One of "total", "background", or "triggering".
#' @param cell_size Optional positive tile width and height. If NULL,
#'   each dimension is inferred from the corresponding coordinate spacing.
#'   Supply this argument when a coordinate axis has only one location.
#'
#' @return A ggplot object.
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
#'                                 days = c(100, 150))
#' plot_intensity_snapshot(grid, time = 150, component = "total")
#' }
#' @export
plot_intensity_snapshot <- function(
    grid,
    time,
    component = c("total", "background", "triggering"),
    cell_size = NULL) {

  component <- match.arg(component)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  required <- c("days", "x", "y", "inside", "snapshots")

  if (!is.list(grid) || !all(required %in% names(grid))) {
    stop(
      "'grid' must be an output of background_trigger_grid().",
      call. = FALSE
    )
  }

  if (!is.numeric(time) || length(time) != 1L || !is.finite(time)) {
    stop("'time' must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(grid$days) || any(!is.finite(grid$days)) ||
      anyDuplicated(grid$days) ||
      length(grid$snapshots) != length(grid$days)) {
    stop("The grid contains invalid times or snapshots.", call. = FALSE)
  }

  snapshot_index <- match(time, grid$days)

  if (is.na(snapshot_index)) {
    stop(
      "'time' is not available in 'grid$days'.",
      call. = FALSE
    )
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

  infer_spacing <- function(z) {
    if (length(z) < 2L) {
      stop(
        "Supply 'cell_size' when a grid axis has only one location.",
        call. = FALSE
      )
    }

    spacing <- diff(z)
    tolerance <- 1e-7 * max(1, abs(spacing[1L]))

    if (any(abs(spacing - spacing[1L]) > tolerance)) {
      stop("Grid coordinates must be regularly spaced.", call. = FALSE)
    }

    spacing[1L]
  }

  if (is.null(cell_size)) {
    tile_width <- infer_spacing(grid$x)
    tile_height <- infer_spacing(grid$y)
  } else {
    if (!is.numeric(cell_size) || length(cell_size) != 1L ||
        !is.finite(cell_size) || cell_size <= 0) {
      stop("'cell_size' must be a positive finite number.",
           call. = FALSE)
    }

    tile_width <- cell_size
    tile_height <- cell_size
  }

  snapshot <- grid$snapshots[[snapshot_index]]
  expected_dim <- c(length(grid$x), length(grid$y))

  extract_matrix <- function(name) {
    z <- snapshot[[name]]

    if (!is.matrix(z) || !is.numeric(z) ||
        !identical(dim(z), expected_dim)) {
      stop("The snapshot contains invalid matrix dimensions or values.",
           call. = FALSE)
    }

    z
  }

  intensity <- switch(
    component,
    background = extract_matrix("bg"),
    triggering = extract_matrix("tr"),
    total = extract_matrix("bg") + extract_matrix("tr")
  )

  coordinates <- expand.grid(x = grid$x, y = grid$y)

  if (!is.logical(grid$inside) ||
      length(grid$inside) != nrow(coordinates) ||
      anyNA(grid$inside)) {
    stop("The grid contains an invalid inside-region mask.",
         call. = FALSE)
  }

  if (!any(grid$inside)) {
    stop("No grid locations are inside the study region.",
         call. = FALSE)
  }

  coordinates$intensity <- as.vector(intensity)
  plot_data <- coordinates[grid$inside, , drop = FALSE]

  if (any(!is.finite(plot_data$intensity)) ||
      any(plot_data$intensity < 0)) {
    stop("Intensity values must be finite and non-negative.",
         call. = FALSE)
  }

  upper_limit <- max(plot_data$intensity)

  if (upper_limit == 0) {
    upper_limit <- 1
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      fill = .data$intensity
    )
  ) +
    ggplot2::geom_tile(
      width = tile_width,
      height = tile_height
    ) +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, upper_limit),
      name = "Intensity"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste(
        tools::toTitleCase(component),
        "intensity"
      ),
      subtitle = paste(
        "Model time:",
        format(time, trim = TRUE),
        "days"
      ),
      x = "x",
      y = "y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

#' @importFrom rlang .data
NULL

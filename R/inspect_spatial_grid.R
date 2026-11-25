#' Inspect a candidate spatial grid
#'
#' Construct and optionally plot a regular grid over a spatial window
#' to inspect a candidate cell size before fitting a model.
#'
#'
#' @param region An \code{sf} or \code{sfc} polygonal spatial window
#'   in a Cartesian coordinate system with units in metres.
#' @param cell_size A positive number giving the grid cell side
#'   length in metres.
#' @param points Optional \code{sf} or \code{sfc} point geometries
#'   to overlay. Must use the same coordinate reference system
#'   as \code{region}.
#' @param plot Logical. Whether to plot the window, grid and points.
#'   Default \code{TRUE}.
#' @param verbose Logical. Whether to print the window dimensions
#'   and number of retained cells. Default \code{TRUE}.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{grid}{An \code{sf} object containing the retained grid cells.}
#'   \item{cell_size}{The requested cell side length in metres.}
#'   \item{width_km, height_km}{
#'     Width and height of the window's bounding box in kilometres.
#'   }
#'   \item{n_cells}{Number of retained grid cells.}
#' }
#'
#' @details
#' This function is a descriptive aid, not an automatic cell-size
#' selector.
#'
#' @examples
#' data("milwaukee")
#' data("milwaukee_county")
#'
#' inspection <- inspect_spatial_grid(
#'   region = milwaukee_county,
#'   cell_size = 500,
#'   points = milwaukee
#' )
#'
#' inspection$n_cells
#'
#' @export
inspect_spatial_grid <- function(region, cell_size, points = NULL,
                                 plot = TRUE, verbose = TRUE) {

  # Validate inputs
  if (!inherits(region, "sf") && !inherits(region, "sfc")) {
    stop("`region` must be an sf or sfc object.", call. = FALSE)
  }

  if (!is.numeric(cell_size) ||
      length(cell_size) != 1L ||
      !is.finite(cell_size) ||
      cell_size <= 0) {
    stop("`cell_size` must be a single positive number.", call. = FALSE)
  }

  region_geom <- sf::st_geometry(region)

  if (length(region_geom) == 0L ||
      any(sf::st_is_empty(region_geom))) {
    stop("`region` must contain non-empty geometries.", call. = FALSE)
  }

  if (!all(sf::st_geometry_type(region_geom) %in%
           c("POLYGON", "MULTIPOLYGON"))) {
    stop("`region` must contain polygon geometries.", call. = FALSE)
  }

  crs <- sf::st_crs(region_geom)

  if (is.na(crs) || !isFALSE(sf::st_is_longlat(region_geom))) {
    stop(
      "`region` must have a Cartesian CRS with units in metres.",
      call. = FALSE
    )
  }

  crs_units <- crs$units_gdal

  if (is.null(crs_units) ||
      length(crs_units) != 1L ||
      is.na(crs_units) ||
      !tolower(crs_units) %in% c("metre", "meter", "metres", "meters", "m")) {
    stop("`region` must use metres as its coordinate units.", call. = FALSE)
  }

  if (!all(sf::st_is_valid(region_geom) %in% TRUE)) {
    stop(
      "`region` contains invalid geometries; repair them first.",
      call. = FALSE
    )
  }

  if (!is.null(points)) {
    if (!inherits(points, "sf") && !inherits(points, "sfc")) {
      stop("`points` must be an sf or sfc object.", call. = FALSE)
    }

    if (!isTRUE(sf::st_crs(points) == crs)) {
      stop(
        "`points` and `region` must use the same CRS.",
        call. = FALSE
      )
    }

    if (!all(sf::st_geometry_type(points) == "POINT")) {
      stop("`points` must contain POINT geometries.", call. = FALSE)
    }
  }

  # Window dimensions
  bb <- sf::st_bbox(region_geom)

  width_km <- unname(bb["xmax"] - bb["xmin"]) / 1000
  height_km <- unname(bb["ymax"] - bb["ymin"]) / 1000

  # Candidate grid
  grid_geom <- sf::st_make_grid(
    region_geom,
    cellsize = cell_size,
    what = "polygons"
  )

  keep <- lengths(
    sf::st_intersects(grid_geom, region_geom)
  ) > 0L

  grid <- sf::st_sf(
    cell_id = seq_len(sum(keep)),
    geometry = grid_geom[keep]
  )

  # Report
  if (isTRUE(verbose)) {
    message(sprintf(
      paste0(
        "Window bounding box: %.2f x %.2f km\n",
        "Cell size: %g m\n",
        "Intersecting cells: %d"
      ),
      width_km, height_km, cell_size, nrow(grid)
    ))
  }

  # Plot
  if (isTRUE(plot)) {
    graphics::plot(
      region_geom,
      col = "grey95",
      border = "grey40",
      main = sprintf("Candidate grid: %g m cells", cell_size)
    )

    graphics::plot(
      sf::st_geometry(grid),
      add = TRUE,
      border = "grey75"
    )

    if (!is.null(points) && length(sf::st_geometry(points)) > 0L) {
      graphics::plot(
        sf::st_geometry(points),
        add = TRUE,
        pch = 16,
        cex = 0.2
      )
    }
  }

  invisible(list(
    grid = grid,
    cell_size = cell_size,
    width_km = width_km,
    height_km = height_km,
    n_cells = nrow(grid)
  ))
}

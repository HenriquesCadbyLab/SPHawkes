#' Convert a semihawkes background grid to an sf polygon grid
#'
#' Turns the rectangular `background$x`/`background$y`/`background$basevalue`
#' grid produced by `background_spatial()` into an `sf` object with one
#' square POLYGON per grid cell, so it can be cropped against `region` with
#' `sf::st_intersection()` and plotted with `ggplot2::geom_sf()`.
#'
#' @param background The `background` component of a `semihawkes_fit`
#'   object (has `$x`, `$y`, `$basevalue`, `$marks`).
#' @param region An `sf`/`sfc` POLYGON/MULTIPOLYGON, used only to supply
#'   the CRS.
#' @return An `sf` object with columns `x`, `y`, `basevalue`, and one square
#'   POLYGON geometry per grid cell.
#' @keywords internal
background_to_grid_sf <- function(background, region) {
  x_vals <- background$x
  y_vals <- background$y

  cellsize_x <- min(diff(sort(unique(x_vals))))
  cellsize_y <- min(diff(sort(unique(y_vals))))

  vals <- background$basevalue
  vals[background$marks <= 0] <- NA_real_

  grid_df <- expand.grid(x = x_vals, y = y_vals)
  grid_df$basevalue <- as.vector(vals)

  make_cell <- function(x, y) {
    sf::st_polygon(list(rbind(
      c(x - cellsize_x / 2, y - cellsize_y / 2),
      c(x + cellsize_x / 2, y - cellsize_y / 2),
      c(x + cellsize_x / 2, y + cellsize_y / 2),
      c(x - cellsize_x / 2, y + cellsize_y / 2),
      c(x - cellsize_x / 2, y - cellsize_y / 2)
    )))
  }

  geoms <- mapply(make_cell, grid_df$x, grid_df$y, SIMPLIFY = FALSE)
  sf::st_sf(
    grid_df,
    geometry = sf::st_sfc(
      geoms,
      crs = sf::st_crs(region)
    ),
    agr = c(
      x = "constant",
      y = "constant",
      basevalue = "constant"
    )
  )
}

#' Crop a semihawkes spatial background to its study-area sf polygon
#'
#' Restricts the fitted spatial-background grid to the `region` the model
#' was fitted on, instead of the full rectangular estimation grid, by
#' intersecting each grid cell (as a small square polygon) with `region`.
#'
#' @details
#' `sf::st_intersection()` between polygons can occasionally return mixed
#' geometry types (e.g. a `POLYGON` clipped down to a `LINESTRING` or
#' `POINT` when it only grazes the boundary, or a `GEOMETRYCOLLECTION`).
#' `ggplot2::geom_sf()`/`fortify()` cannot handle a geometry column with
#' mixed/non-standard types, so we explicitly keep only (MULTI)POLYGON
#' results and cast them to POLYGON/MULTIPOLYGON.
#'
#' @param background The `background` component of a `semihawkes_fit`
#'   object.
#' @param region An `sf`/`sfc` POLYGON/MULTIPOLYGON study-area boundary,
#'   typically `fit$events$region`.
#' @param crop Currently unused; kept for interface compatibility. Cells
#'   outside `region` are always dropped (they have no intersection).
#' @return An `sf` object with columns `x`, `y`, `basevalue`, and one
#'   (possibly clipped) polygon geometry per grid cell that overlaps
#'   `region`, ready for `ggplot2::geom_sf(aes(fill = basevalue))`.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' bg <- background_spatial(ev, cell_size = 1000)
#' cells <- crop_background_spatial(bg, ev$region)
#' head(cells)
#' @export
crop_background_spatial <- function(background, region, crop = TRUE) {
  grid_sf <- background_to_grid_sf(background, region)
  region_valid <- sf::st_make_valid(region)
  region_union <- sf::st_sf(geometry = sf::st_union(region_valid))

  clipped <- sf::st_intersection(grid_sf, region_union)

  geom_types <- sf::st_geometry_type(clipped)
  clipped <- clipped[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
  clipped <- sf::st_cast(clipped, "MULTIPOLYGON")

  clipped
}

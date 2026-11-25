#' @keywords internal
#' @noRd
.check_sf_points <- function(points, arg = "points") {
  if (!inherits(points, "sf")) {
    stop(sprintf("`%s` must be an object of class 'sf' (got '%s').", arg,
                 paste(class(points), collapse = "/")), call. = FALSE)
  }
  geom_types <- unique(as.character(sf::st_geometry_type(points)))
  if (!all(geom_types %in% c("POINT"))) {
    stop(sprintf("`%s` must have POINT geometry only (found: %s).", arg,
                 paste(geom_types, collapse = ", ")), call. = FALSE)
  }
  if (is.na(sf::st_crs(points))) {
    stop(sprintf("`%s` has no coordinate reference system (CRS) set. ",
                 arg), " Set one with sf::st_set_crs() before proceeding.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.check_sf_polygon <- function(region, arg = "region") {
  if (!inherits(region, "sf") && !inherits(region, "sfc")) {
    stop(sprintf("`%s` must be an 'sf' or 'sfc' object (got '%s').", arg,
                 paste(class(region), collapse = "/")), call. = FALSE)
  }
  geom_types <- unique(as.character(sf::st_geometry_type(region)))
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop(sprintf("`%s` must have POLYGON/MULTIPOLYGON geometry (found: %s).",
                 arg, paste(geom_types, collapse = ", ")), call. = FALSE)
  }
  if (is.na(sf::st_crs(region))) {
    stop(sprintf("`%s` has no coordinate reference system (CRS) set.", arg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
# Guard clause: errors out if `x` is in a geographic (lon/lat) CRS instead
# of a projected one. Every kernel bandwidth and distance calculation in
# this package assumes linear units (metres), so geographic coordinates
# would silently produce wrong results rather than an obvious error. Stops
# immediately with a clear message (which CRS was found, and how to fix it
# via sf::st_transform()) instead of failing confusingly downstream.
.check_projected_crs <- function(x, arg = "points") {
  crs <- sf::st_crs(x)
  is_geographic <- isTRUE(sf::st_is_longlat(x))
  if (is_geographic) {
    stop(sprintf(
      paste0(
        "`%s` is in a geographic (longitude/latitude) CRS (%s). ",
        "This package's kernel bandwidths and distances are defined in ",
        "linear units (typically metres), so `%s` must be transformed to a ",
        "projected CRS first, e.g. sf::st_transform(%s, crs = <EPSG code>)."
      ),
      arg, crs$input %||% crs$Name %||% "unknown", arg, arg
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Warn when the working CRS is projected but not in metres
#'
#' `.check_projected_crs()` already rejects geographic (lon/lat) input. This
#' is the softer companion check: a projected CRS in feet, US survey feet or
#' kilometres will work, but every default in the package (bandwidths,
#' `cell_size`, trigger window distances, kernel grid extents) was chosen for
#' metres, so those defaults will be badly wrong in scale. Warn rather than
#' error, since the values are all configurable.
#'
#' @keywords internal
#' @noRd
.check_crs_units <- function(x, arg = "region") {
  units <- tryCatch(sf::st_crs(x)$units, error = function(e) NULL)
  if (is.null(units) || is.na(units)) return(invisible(TRUE))
  if (!identical(units, "m")) {
    warning(sprintf(
      paste0(
        "`%s` uses a projected CRS whose linear unit is '%s', not metres. ",
        "The package's default bandwidths, `cell_size`, trigger windows and ",
        "kernel grid extents are all expressed in the CRS's linear units and ",
        "were chosen for metres, so they will be wrong by a large factor ",
        "here. Either transform to a metre-based CRS, or set every distance ",
        "argument explicitly in '%s'."
      ), arg, units, units
    ), call. = FALSE)
  }
  invisible(TRUE)
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) != 1L || is.na(a)) b else a
}

#' A short, printable label for a CRS
#'
#' `sf::st_crs()$input` is whatever the CRS was built from, which for a CRS
#' constructed from WKT is the entire multi-line WKT string — far too long to
#' drop into a `print()` method. Prefer the SRID (`"EPSG:32717"`), then the
#' human-readable name, then a truncated first line of the input.
#'
#' @keywords internal
#' @noRd
.crs_label <- function(crs) {
  crs <- tryCatch(sf::st_crs(crs), error = function(e) NULL)
  if (is.null(crs) || is.na(crs)) return("unknown")

  short <- function(v) {
    if (is.null(v) || length(v) != 1L || is.na(v)) return(NULL)
    v <- trimws(strsplit(v, "\n", fixed = TRUE)[[1]][1])
    if (!nzchar(v) || nchar(v) > 60) NULL else v
  }

  short(crs$srid) %||% short(crs$Name) %||% short(crs$input) %||% "unknown"
}

#' Regular sequence from `from` that is guaranteed to reach `to`
#'
#' [seq()] with `by` stops at or before `to`, so a grid built with it falls
#' short whenever `by` does not divide `to - from` exactly, leaving a strip of
#' the study area uncovered. This always covers the full extent, overshooting
#' by less than one step when necessary.
#'
#' @param from,to Numeric endpoints, `from < to`.
#' @param by Positive step.
#' @return Numeric vector starting at `from`, spaced by `by`, ending at or
#'   after `to`.
#' @keywords internal
#' @noRd
.cover_seq <- function(from, to, by) {
  # The tolerance stops floating-point drift from adding a spurious extra
  # cell when `by` divides `to - from` exactly.
  n <- max(ceiling((to - from) / by - 1e-9), 1)
  from + by * seq.int(0L, n)
}

#' Area represented by one cell of a background grid
#'
#' Every spatial integral in the package is a Riemann sum over
#' `background_spatial()`'s grid, of the form `sum(values * marks) *
#' cell_area`. Calibrating `cell_area` so that the in-region cells tile the
#' region exactly makes those sums agree with each other and makes the
#' integral of the (mean-1-normalised) background come out as the region's
#' true area, whatever `cell_size` was chosen.
#'
#' This replaces an earlier `bounding box area / number of grid points`,
#' which assumed the grid lined up with the bounding box and that the region
#' filled it.
#'
#' @param spatial A `"semihawkes_spatial_base"` object.
#' @param region The study-region geometry (`events$region`).
#' @return A single positive number, in squared CRS units.
#' @keywords internal
#' @noRd
.grid_cell_area <- function(spatial, region) {
  n_in <- sum(spatial$marks > 0)
  if (n_in == 0) {
    stop("No background grid cell falls inside the study region; `cell_size` ",
         "is too coarse for this region.", call. = FALSE)
  }
  as.numeric(sf::st_area(region)) / n_in
}

#' Divide, treating a zero or non-finite denominator as "no information"
#'
#' The edge corrections in [trigger_reestimate_spatial()] and
#' [trigger_reestimate_temporal()] divide an observed density by how often
#' each offset could have been observed at all. That denominator is zero for
#' offsets no pair of events in the study region could ever realise -- every
#' offset larger than the region, for instance. The numerator is then
#' necessarily zero too, so the honest value of the ratio is zero, but `0 / 0`
#' gives `NaN`, and a single `NaN` entering an FFT smoother spreads across the
#' whole kernel. The result used to be an all-`NA` triggering kernel, an
#' infinite branching ratio, and a fit that failed with `missing value where
#' TRUE/FALSE needed` from the convergence check.
#'
#' @param num,den Numeric vectors or matrices of the same shape.
#' @return `num / den` where `den` is finite and positive, and `0` elsewhere,
#'   with the shape of `num` preserved.
#' @keywords internal
#' @noRd
.ratio_or_zero <- function(num, den) {
  out <- numeric(length(num))
  ok <- is.finite(den) & den > 0 & is.finite(num)
  out[ok] <- num[ok] / den[ok]
  if (!is.null(dim(num))) dim(out) <- dim(num)
  out
}

#' Normalise a re-estimated kernel, failing clearly if it cannot be done
#'
#' @param basevalue Numeric vector or matrix of unnormalised kernel values.
#' @param integral The kernel's integral, from [simpson()] / [simpson.2D()].
#' @param what Name used in the error message.
#' @return `basevalue / integral`.
#' @keywords internal
#' @noRd
.normalise_kernel <- function(basevalue, integral, what) {
  if (!is.finite(integral) || integral <= 0) {
    stop("The re-estimated ", what, " triggering kernel integrates to ",
         format(integral), ", so it cannot be normalised. This usually means ",
         "no event pair carried usable weight at this iteration -- too few ",
         "events, or a triggering kernel whose extent is far larger than the ",
         "study region.", call. = FALSE)
  }
  basevalue / integral
}

#' Save and restore the global random-number state
#'
#' Functions with a `seed` argument seed the RNG for their own draws, then put
#' the caller's RNG state back on exit, so that supplying `seed` makes that
#' call reproducible without changing any random numbers the caller draws
#' afterwards. Use as
#' `old <- .get_rng_state(); on.exit(.restore_rng_state(old), add = TRUE);
#' set.seed(seed)`.
#'
#' @param state The value returned by `.get_rng_state()`.
#' @keywords internal
#' @noRd
.get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
.restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

#' Check that an argument is a single TRUE or FALSE
#'
#' @param x The value to check.
#' @param arg The argument name, for the error message.
#' @keywords internal
#' @noRd
.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Drop a local engineering CRS before plotting with geom_sf()
#'
#' geom_sf()/coord_sf() run every layer through sf::st_transform(), and
#' GDAL/PROJ cannot build any transformation for a local engineering CRS
#' (such as the "Local Cartesian coordinates" of the bundled `crimes` data),
#' not even to itself, so the plot fails with "transformation not available".
#' Such coordinates have no geographic meaning anyway, so they are plotted in
#' plain map units. Any other CRS is kept, so maps of projected data are
#' unchanged.
#'
#' @param x An `sf` or `sfc` object.
#' @return `x`, with its CRS removed if it was an engineering CRS.
#' @keywords internal
#' @noRd
.plot_safe_crs <- function(x) {
  crs <- sf::st_crs(x)
  if (!is.na(crs) && grepl("^\\s*ENGCRS", crs$wkt)) {
    x <- sf::st_set_crs(x, NA)
  }
  x
}

#' Point-in-region test for many points at once
#'
#' Builds a function that says which points intersect `region` (inside it or
#' on its boundary): the same answer as
#' `lengths(sf::st_intersects(points, region)) > 0`, but much faster for the
#' millions of points `trigger_spatial_marks()` tests. Most points are tested
#' with spatstat's `inside.owin()` against two copies of the region, one
#' shrunk and one grown by `eps`. A point inside the shrunk copy is certainly
#' inside the region, and a point outside the grown copy is certainly
#' outside. Only the points in the thin band between them, within `eps` of
#' the boundary, are tested with `sf::st_intersects()` itself, so points on
#' or next to the boundary are classified exactly as before.
#'
#' If buffering cannot produce copies that bracket the region (for example a
#' region narrower than `2 * eps`), every point is tested with
#' `sf::st_intersects()`.
#'
#' @param region An `sf`/`sfc` polygon layer in a projected CRS.
#' @param eps Width of the band, in the CRS's units: far larger than
#'   floating-point error, far smaller than the region's features.
#' @return A function of `x`, `y` (numeric vectors) returning a logical
#'   vector.
#' @keywords internal
#' @noRd
.region_point_test <- function(region, eps) {
  crs <- sf::st_crs(region)
  sf_test <- function(x, y) {
    if (length(x) == 0) return(logical(0))
    pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = crs)
    lengths(sf::st_intersects(pts, region)) > 0
  }

  whole <- sf::st_union(sf::st_geometry(region))
  # Mitred corners keep the copies' outlines as simple as the region's own.
  shrunk <- sf::st_buffer(whole, -eps, joinStyle = "MITRE", mitreLimit = 2)
  grown  <- sf::st_buffer(whole,  eps, joinStyle = "MITRE", mitreLimit = 2)
  brackets <- !any(sf::st_is_empty(shrunk)) &&
    isTRUE(sf::st_covers(whole, shrunk, sparse = FALSE)[1, 1]) &&
    isTRUE(sf::st_covers(grown, whole, sparse = FALSE)[1, 1])
  if (!brackets) return(sf_test)

  w_shrunk <- spatstat.geom::as.owin(shrunk)
  w_grown  <- spatstat.geom::as.owin(grown)
  function(x, y) {
    inside <- spatstat.geom::inside.owin(x, y, w_shrunk)
    band <- which(!inside & spatstat.geom::inside.owin(x, y, w_grown))
    inside[band] <- sf_test(x[band], y[band])
    inside
  }
}

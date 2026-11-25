#' Monte Carlo integral of each event's isotropic Gaussian kernel over the
#' study region
#'
#' For every event, its background kernel (an isotropic Gaussian centred on
#' the event with standard deviation `bandwidth[i]`) generally has some of
#' its mass falling outside the study region boundary. This computes the
#' fraction of that mass that falls *inside* `region`, by Monte Carlo
#' integration over `n_mc` random points drawn uniformly in the region. It
#' is used to edge-correct the background-density kernel estimate.
#'
#' @param coorx,coory Numeric vectors of projected event coordinates.
#' @param bandwidth Numeric vector (same length), the per-event standard
#'   deviation from e.g. [adaptive_bandwidth()].
#' @param region An `sf`/`sfc` polygon (already in the same projected CRS as
#'   the coordinates).
#' @param n_mc Number of Monte Carlo points drawn in `region`; default
#'   `3000`, matching the original script.
#' @param cap Values above `cap` are truncated to `1` (default `0.99`), since
#'   a fully-interior kernel should integrate to ~1 and larger Monte Carlo
#'   noise is expected near that ceiling.
#' @param seed Optional random seed, making the Monte Carlo draw
#'   reproducible.
#' @param restore_rng Only used when `seed` is supplied. If `TRUE` (the
#'   default), the global random-number state is put back as it was before
#'   the call, so the seed affects only this call. If `FALSE`, the call leaves
#'   the generator in the state `set.seed(seed)` and its own draws produce,
#'   as earlier versions of the package did -- useful for reproducing results
#'   from scripts that relied on that.
#'
#' @return A numeric vector (same length as `coorx`) with each event's
#'   in-region kernel mass, in `(0, 1]`.
#' @examples
#' xy <- sf::st_coordinates(sim_events)[1:5, ]
#' # Share of a 300 m Gaussian kernel around each event that lies in the region
#' bg_integral(xy[, 1], xy[, 2], bandwidth = rep(300, 5), region = sim_region,
#'             n_mc = 500, seed = 1)
#' @export
bg_integral <- function(coorx, coory, bandwidth, region, n_mc = 3000,
                        cap = 0.99, seed = NULL, restore_rng = TRUE) {
  n <- length(coorx)
  if (length(coory) != n || length(bandwidth) != n) {
    stop("`coorx`, `coory` and `bandwidth` must have the same length.",
         call. = FALSE)
  }
  .check_sf_polygon(region, "region")

  .check_flag(restore_rng, "restore_rng")
  if (!is.null(seed)) {
    if (restore_rng) {
      old_rng <- .get_rng_state()
      on.exit(.restore_rng_state(old_rng), add = TRUE)
    }
    set.seed(seed)
  }

  w <- spatstat.geom::as.owin(region)
  pts <- spatstat.random::rpoint(n_mc, win = w)
  mc_coords <- cbind(pts$x, pts$y)
  area_w <- spatstat.geom::area.owin(w)

  gaussian_density <- function(coords, mean, sigma_diag) {
    diffs <- coords
    diffs[, 1] <- diffs[, 1] - mean[1]
    diffs[, 2] <- diffs[, 2] - mean[2]
    exp_term <- rowSums(diffs^2) / sigma_diag
    (1 / (2 * pi * sigma_diag)) * exp(-0.5 * exp_term)
  }

  out <- vapply(seq_len(n), function(i) {
    mean_i <- c(coorx[i], coory[i])
    sigma_diag <- bandwidth[i]^2
    mean(gaussian_density(mc_coords, mean = mean_i, sigma_diag = sigma_diag)) * area_w
  }, numeric(1))

  ifelse(out > cap, 1, out)
}

#' Add the background-kernel edge-correction integral to a
#' `semihawkes_events` object
#'
#' Convenience wrapper around [bg_integral()]. Requires `x$events$bandwidth`
#' to already exist (see [add_adaptive_bandwidth()]), and additionally drops
#' any event that falls outside `x$region` after re-checking containment.
#'
#' @param x A `semihawkes_events` object with a `bandwidth` column already
#'   added.
#' @param ... Passed on to [bg_integral()] (`n_mc`, `cap`, `seed`,
#'   `restore_rng`).
#' @return `x`, with a `bg.integral` column added to `x$events`.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' summary(ev$events$bg.integral)
#' @export
add_bg_integral <- function(x, ...) {
  if (!inherits(x, "semihawkes_events")) {
    stop("`x` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  if (is.null(x$events$bandwidth)) {
    stop("`x$events$bandwidth` is missing; call `add_adaptive_bandwidth()` first.",
         call. = FALSE)
  }
  x$events$bg.integral <- bg_integral(x$events$coorx, x$events$coory,
                                      x$events$bandwidth, x$region, ...)

  w <- spatstat.geom::as.owin(x$region)
  inside <- spatstat.geom::inside.owin(x = x$events$coorx, y = x$events$coory, w = w)
  if (any(!inside)) {
    x$events <- x$events[inside, ]
    x$n_kept <- nrow(x$events)
  }
  x
}

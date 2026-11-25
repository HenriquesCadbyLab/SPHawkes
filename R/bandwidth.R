#' Adaptive (nearest-neighbour) spatial bandwidth for each event
#'
#' For every point, the bandwidth is the distance to its `k`-th nearest
#' neighbour among neighbours at least `min_dist` away (this skips
#' near-duplicate/co-located points, e.g. events geocoded to the same
#' address), floored at `min_bandwidth`.
#'
#' @param x,y Numeric vectors of projected event coordinates (same length).
#' @param k Which nearest neighbour (beyond `min_dist`) to use; default `5`.
#' @param min_dist Minimum distance for a neighbour to be considered
#'   (neighbours closer than this are ignored), in the linear units of the working CRS (metres for a UTM projection);
#'   default `150`, i.e. 150 m.
#' @param min_bandwidth Floor applied to every bandwidth, in the linear units of the working CRS (metres for a UTM projection); default `250`, i.e. 250 m.
#' @param k_search How many nearest neighbours to search over before
#'   filtering by `min_dist`; increase this if a warning about insufficient
#'   neighbours is raised. Default `50`.
#'
#' @return A numeric vector of bandwidths, one per input point.
#' @examples
#' xy <- sf::st_coordinates(sim_events)
#' bw <- adaptive_bandwidth(xy[, 1], xy[, 2], k = 5, min_dist = 150,
#'                          min_bandwidth = 300)
#' summary(bw)
#' @export
adaptive_bandwidth <- function(x, y, k = 5, min_dist = 150,
                               min_bandwidth = 250, k_search = 50) {
  n <- length(x)
  if (length(y) != n) stop("`x` and `y` must have the same length.", call. = FALSE)
  if (n <= k) {
    warning("Fewer points than `k`; returning `min_bandwidth` for all points.",
            call. = FALSE)
    return(rep(min_bandwidth, n))
  }

  k_search <- min(k_search, n - 1)

  win <- spatstat.geom::owin(range(x) + c(-1, 1), range(y) + c(-1, 1))
  X <- spatstat.geom::ppp(x, y, window = win, check = FALSE)

  # nndist(..., k = 1:k_search) returns an n x k_search matrix, each row
  # already sorted ascending (1st NN, 2nd NN, ...).
  dmat <- spatstat.geom::nndist(X, k = seq_len(k_search))
  if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = 1)

  insufficient <- FALSE
  bandwidth <- vapply(seq_len(n), function(i) {
    d <- dmat[i, ]
    d <- d[d >= min_dist]
    if (length(d) < k) {
      # Fall back to an exact search for this point only.
      dd <- sqrt((x - x[i])^2 + (y - y[i])^2)[-i]
      dd <- sort(dd[dd >= min_dist])
      if (length(dd) < k) {
        insufficient <<- TRUE
        return(NA_real_)
      }
      dd[k]
    } else {
      d[k]
    }
  }, numeric(1))

  if (insufficient) {
    warning(sprintf(
      "Fewer than %d neighbours found beyond `min_dist` = %s for some points; ",
      k, min_dist), "using `min_bandwidth` for those points. Consider lowering ",
      "`min_dist`/`k` or increasing `k_search`.", call. = FALSE)
    bandwidth[is.na(bandwidth)] <- min_bandwidth
  }

  pmax(min_bandwidth, bandwidth)
}

#' Add an adaptive bandwidth column to a `semihawkes_events` object
#'
#' Convenience wrapper around [adaptive_bandwidth()] that operates directly
#' on the `events` data.frame produced by [prepare_events()].
#'
#' @param x A `semihawkes_events` object (from [prepare_events()]).
#' @param ... Passed on to [adaptive_bandwidth()] (`k`, `min_dist`,
#'   `min_bandwidth`, `k_search`).
#' @return `x`, with a `bandwidth` column added to `x$events`.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' summary(ev$events$bandwidth)
#' @export
add_adaptive_bandwidth <- function(x, ...) {
  if (!inherits(x, "semihawkes_events")) {
    stop("`x` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  x$events$bandwidth <- adaptive_bandwidth(x$events$coorx, x$events$coory, ...)
  x
}

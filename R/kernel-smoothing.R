# Numerical-integration and kernel-smoothing primitives used throughout the
# EM algorithm: Simpson's rule integration, weighted histograms, and Gaussian
# kernel smoothing via direct convolution or FFT.
#
# Adapted from the supplementary R code accompanying D'Angelo et al. (2024).
# Modifications in spHawkes include namespace qualification, documentation,
# histogram boundary handling, and corrections for rectangular grids.
#
# @reference
# D'Angelo, N., Payares, D., Adelfio, G., and Mateu, J. (2024).
# Self-exciting point process modelling of crimes on linear networks.
# Statistical Modelling, 24(2), 139-168. \doi{10.1177/1471082X221094146}

#' Simpson's rule integral of an equally-spaced function
#'
#' Approximates \eqn{\int f(x)\,dx} from function values `value` sampled on
#' an equally-spaced grid with step `dx`, using (composite) Simpson's rule.
#'
#' @param value Numeric vector of function values on an equally-spaced grid.
#'   Simpson's rule needs an odd number of points (an even number of
#'   intervals), at least 3; any other length is an error.
#' @param dx Grid spacing; default `1`.
#'
#' @return A single numeric value, the approximate integral.
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#' @examples
#' # The standard normal density integrates to (almost exactly) 1
#' x <- seq(-4, 4, by = 0.1)
#' simpson(stats::dnorm(x), dx = 0.1)
#' @export
simpson <- function(value, dx = 1) {
  m <- length(value)
  # The weights 1, 4, 2, 4, ..., 2, 4, 1 only fit an odd number of points.
  # For an even number, `rep()` truncates and the shorter weight vector is
  # recycled, giving a wrong integral with only a recycling warning.
  if (m < 3 || m %% 2 == 0) {
    stop("Simpson's rule needs an odd number of equally spaced points, at ",
         "least 3; `value` has ", m, ".", call. = FALSE)
  }
  bx <- c(1, rep(c(4, 2), (m - 3) / 2), 4, 1)
  sum(value * bx) * dx / 3.0
}

#' Simpson's rule double integral of an equally-spaced function
#'
#' 2D analogue of [simpson()]: approximates \eqn{\int\int f(x,y)\,dx\,dy}
#' from a matrix of function values sampled on an equally-spaced grid in
#' both dimensions.
#'
#' @param value Numeric matrix of function values, rows indexing the first
#'   dimension (step `dx`) and columns the second (step `dy`). As for
#'   [simpson()], both the number of rows and the number of columns must be
#'   odd and at least 3.
#' @param dx,dy Grid spacing along rows/columns; default `1` each.
#'
#' @return A single numeric value, the approximate double integral.
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#' @examples
#' g <- seq(-4, 4, by = 0.1)
#' z <- outer(stats::dnorm(g), stats::dnorm(g))
#' simpson.2D(z, dx = 0.1, dy = 0.1)
#' @export
simpson.2D <- function(value, dx = 1, dy = 1) {
  m <- nrow(value)
  n <- ncol(value)
  # See simpson(): the weights only fit an odd number of points. Here an even
  # number makes the weight vectors one short, so the matrix product fails.
  if (m < 3 || n < 3 || m %% 2 == 0 || n %% 2 == 0) {
    stop("Simpson's rule needs an odd number of equally spaced points, at ",
         "least 3, in each direction; `value` is ", m, " x ", n, ".",
         call. = FALSE)
  }

  bx <- c(1, rep(c(4, 2), (m - 3) / 2), 4, 1)
  by <- c(1, rep(c(4, 2), (n - 3) / 2), 4, 1)

  sum((bx %*% value %*% by) * dx * dy / 9.0)
}

# Checks a kernel grid's span and step for trigger_temporal_init() and
# trigger_spatial_init(): `span / step` must be a whole number -- and, with
# `even = TRUE`, an even one -- so that the grid ends exactly at `span` and has
# the odd number of points simpson() and simpson.2D() need. `grid` completes
# the error message's description of the grid. Returns the number of steps.
.grid_steps <- function(span, step, span_arg, step_arg, grid, even = FALSE) {
  positive <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v > 0
  if (!positive(span) || !positive(step)) {
    stop(sprintf("`%s` and `%s` must each be a single positive number.",
                 span_arg, step_arg), call. = FALSE)
  }
  k <- span / step
  whole <- abs(k - round(k)) <= 1e-8 * max(1, k)
  if (!whole || (even && round(k) %% 2 != 0)) {
    stop(sprintf(paste0(
      "`%s / %s` must be %s whole number, so that the grid %s and has the ",
      "odd number of points Simpson's rule needs. `%s = %s` and `%s = %s` ",
      "give %s."),
      span_arg, step_arg, if (even) "an even" else "a", grid,
      span_arg, format(span), step_arg, format(step), format(signif(k, 6))),
      call. = FALSE)
  }
  invisible(round(k))
}

#' Weighted 1D histogram
#'
#' Bins `x` into `breaks` and sums `weight` (instead of counting) within
#' each bin, returning both counts and a weighted density (counts divided
#' by total weight and by bin width) plus bin midpoints — the building
#' block for the non-parametric weekly/monthly/trend estimators.
#'
#' @param x Numeric vector of values to bin.
#' @param weight Numeric vector of weights, same length as `x` (e.g. the
#'   EM background/trigger probabilities).
#' @param breaks Numeric vector of bin edges; must span the range of `x`
#'   (an error otherwise).
#'
#' @return A list with `breaks`, `counts` (summed weight per bin),
#'   `density` (`counts` normalised to integrate to 1 over `breaks`), and
#'   `mids` (bin midpoints).
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#'
#' @section A note on the name:
#' Despite the dots, this is an ordinary function, **not** an S3 method for
#' the [graphics::hist()] generic — the name is inherited from the original
#' D'Angelo et al. (2024) script. It has to be exported with an explicit
#' `@rawNamespace` directive because roxygen2 otherwise reads `hist.weighted`
#' as `hist` dispatching on class `"weighted"` and registers it as a method,
#' which both hides the function from users and makes `hist()` misbehave on
#' any object that happens to carry that class.
#'
#' @examples
#' x <- c(0.2, 0.4, 1.3, 2.7, 2.9)
#' w <- c(1, 2, 1, 1, 3)
#' h <- hist.weighted(x, w, breaks = 0:3)
#' h$counts
#' h$density
#' @rawNamespace export(hist.weighted)
hist.weighted <- function(x, weight, breaks) {
  # The original wrapped this check in `try()`, which printed the error and
  # carried on, so out-of-range values surfaced later as an unrelated length
  # mismatch. Stop here instead.
  if (any(x < min(breaks) | x > max(breaks), na.rm = TRUE)) {
    stop(sprintf("`breaks` (%s to %s) do not span the range of `x` (%s to %s).",
                 format(min(breaks)), format(max(breaks)),
                 format(min(x, na.rm = TRUE)), format(max(x, na.rm = TRUE))),
         call. = FALSE)
  }
  nb <- length(breaks)

  # f(break_i) = 0:n_b[i], f(x) = c for x in (break_i, break_[i+1])
  tempfun <- stats::stepfun(breaks, 0:nb)
  mids <- (breaks[2:nb] + breaks[1:(nb - 1)]) * .5
  x <- c(x, mids)
  weight <- c(weight, rep(0, nb - 1))
  z <- tempfun(x)
  # Snap only exact boundary matches into the nearest closed bin (fix, see
  # @section above). Out-of-range `x` has already stopped above: unlike a
  # legitimate boundary point, there is no single correct bin for it.
  z[z == 0 & x == min(breaks)] <- 1L
  z[z == nb & x == max(breaks)] <- nb - 1L
  freq <- tapply(weight, z, sum)
  list(breaks = breaks, counts = freq,
       density = freq / sum(freq) / (breaks[2:nb] - breaks[1:(nb - 1)]),
       mids = mids)
}

#' @rdname hist.weighted
#' @rawNamespace export(hist.weighted.1D)
hist.weighted.1D <- hist.weighted

#' Weighted 2D histogram
#'
#' 2D analogue of [hist.weighted()]: bins `(x, y)` pairs into a grid defined
#' by `x.breaks`/`y.breaks` and sums `weight` within each cell.
#'
#' @param x,y Numeric vectors of coordinates to bin (same length).
#' @param weight Numeric vector of weights, same length as `x`/`y`.
#' @param x.breaks,y.breaks Numeric vectors of bin edges for each dimension;
#'   must span the range of `x`/`y` respectively (an error otherwise).
#'
#' @return A list with `x.breaks`, `y.breaks`, `x.mids`, `y.mids`, `counts`
#'   (matrix of summed weight per cell), and `density` (`counts` normalised
#'   to integrate to 1 over the grid).
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#'
#' @section A note on the name:
#' As with [hist.weighted()], this is an ordinary function rather than an S3
#' method for [graphics::hist()], and needs an explicit `@rawNamespace`
#' export directive for that reason.
#'
#' @examples
#' h <- hist.weighted.2D(
#'   x = c(0.2, 0.7, 1.5), y = c(0.1, 0.9, 1.2), weight = c(1, 2, 1),
#'   x.breaks = 0:2, y.breaks = 0:2
#' )
#' h$density
#' @rawNamespace export(hist.weighted.2D)
hist.weighted.2D <- function(x, y, weight, x.breaks, y.breaks) {
  # As in hist.weighted(): stop rather than report from inside `try()`.
  if (any(x < min(x.breaks) | x > max(x.breaks), na.rm = TRUE) ||
      any(y < min(y.breaks) | y > max(y.breaks), na.rm = TRUE)) {
    stop(sprintf(paste0("`x.breaks` (%s to %s) and `y.breaks` (%s to %s) do ",
                        "not span the range of `x` (%s to %s) and `y` (%s to %s)."),
                 format(min(x.breaks)), format(max(x.breaks)),
                 format(min(y.breaks)), format(max(y.breaks)),
                 format(min(x, na.rm = TRUE)), format(max(x, na.rm = TRUE)),
                 format(min(y, na.rm = TRUE)), format(max(y, na.rm = TRUE))),
         call. = FALSE)
  }
  nbx <- length(x.breaks)
  nby <- length(y.breaks)

  mids.x <- (x.breaks[2:nbx] + x.breaks[1:(nbx - 1)]) * .5
  mids.y <- (y.breaks[2:nby] + y.breaks[1:(nby - 1)]) * .5
  mids.mat.x <- mids.x %o% rep(1, nby - 1)
  mids.mat.y <- rep(1, nbx - 1) %o% mids.y # nolint: fixed, see @section above (was `nby - 1`)
  xx <- c(x, mids.mat.x)
  yy <- c(y, mids.mat.y)
  weight <- c(weight, rep(0, (nbx - 1) * (nby - 1)))

  tempfun.x <- stats::stepfun(x.breaks, 0:nbx)
  tempfun.y <- stats::stepfun(y.breaks, 0:nby)
  tx <- tempfun.x(xx)
  ty <- tempfun.y(yy)
  # Snap only exact boundary matches (same rationale as hist.weighted()'s
  # fix, above); out-of-range points have already stopped.
  tx[tx == 0 & xx == min(x.breaks)] <- 1L
  tx[tx == nbx & xx == max(x.breaks)] <- nbx - 1L
  ty[ty == 0 & yy == min(y.breaks)] <- 1L
  ty[ty == nby & yy == max(y.breaks)] <- nby - 1L
  z <- tx + ty * nbx

  freq <- matrix(tapply(weight, z, sum), nrow = nbx - 1, ncol = nby - 1)

  x1.mat <- x.breaks %o% rep(1, nby)
  y1.mat <- rep(1, nbx) %o% y.breaks # nolint: fixed, see @section above (was `nby`)

  list(x.breaks = x.breaks, y.breaks = y.breaks, x.mids = mids.x, y.mids = mids.y,
       counts = freq,
       density = freq / sum(freq) /
         (x1.mat[2:nbx, 2:nby] - x1.mat[1:(nbx - 1), 1:(nby - 1)]) /
         (y1.mat[2:nbx, 2:nby] - y1.mat[1:(nbx - 1), 1:(nby - 1)]))
}

#' 1D Gaussian kernel smoothing by direct convolution, with reflected boundaries
#'
#' Smooths values `y` sampled at equally-spaced points `x` with a Gaussian
#' kernel of standard deviation `bandwidth`, via [stats::convolve()],
#' extending the series by reflection (`c(y, rev(y))`) at both ends as in the
#' original implementation. See [ker.smooth.fft()] for the periodic-boundary
#' counterpart and why the two are not interchangeable.
#'
#' @param x Numeric vector of equally-spaced grid points.
#' @param y Numeric vector of values at `x` (same length).
#' @param bandwidth Standard deviation of the Gaussian smoothing kernel.
#'
#' @return A numeric vector of length `length(x) + 1`, the smoothed values.
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#' @examples
#' # Smooth a weighted histogram; the result has one value per break
#' breaks <- seq(0, 10, by = 0.1)
#' h <- hist.weighted(c(2, 2.5, 6, 6.2, 6.4), rep(1, 5), breaks = breaks)
#' smooth <- ker.smooth.conv(h$mids, h$density, bandwidth = 0.5)
#' plot(breaks, smooth, type = "l", xlab = "x", ylab = "Smoothed density")
#' @export
ker.smooth.conv <- function(x, y, bandwidth) {
  yy <- c(y, rev(y))

  n <- length(x)
  dx <- x[2] - x[1]
  xx <- ((1:n) * dx - dx / 2)

  ker.val <- stats::dnorm(xx, 0, sd = bandwidth)
  ker.val <- c(ker.val, rev(ker.val))

  smoothed <- stats::convolve(yy, ker.val) * dx

  smoothed[1:(n + 1)]
}

#' 1D Gaussian kernel smoothing by FFT, with periodic boundaries
#'
#' Smooths values `z` sampled at equally-spaced points `x` with a Gaussian
#' kernel of standard deviation `bandwidth`, using an FFT.
#'
#' @section Not interchangeable with [ker.smooth.conv()]:
#' This is **not** simply a faster version of [ker.smooth.conv()]: the two
#' treat the boundaries differently, and generally return different values.
#'
#' * [ker.smooth.conv()] extends the series by *reflection* (`c(y, rev(y))`),
#'   which is the right assumption for a curve that is not periodic — such as
#'   the temporal triggering kernel, whose lag axis simply stops at the
#'   cutoff.
#' * `ker.smooth.fft()` extends the series *periodically* (`c(z, z)`), which
#'   is the right assumption for a genuinely cyclical curve — such as the
#'   day-of-week and month-of-year background components, where the value at
#'   Saturday midnight really does wrap round to Sunday.
#'
#' The package uses each accordingly: [temporal_base_weekly()] and
#' [temporal_base_monthly()] call this function, while
#' [trigger_reestimate_temporal()] calls [ker.smooth.conv()], and
#' [trigger_reestimate_spatial()] calls [ker.smooth.2D.fft()] (which, despite
#' also being FFT-based, reflects rather than wraps). Substituting one for
#' another "for speed" will silently change the estimates near the
#' boundaries.
#'
#' @param x Numeric vector of equally-spaced grid points.
#' @param z Numeric vector of values at `x` (same length).
#' @param bandwidth Standard deviation of the Gaussian smoothing kernel.
#'
#' @return A numeric vector of length `length(z) + 1`: the smoothed values
#'   at the points halfway between successive `x`, and half a grid step
#'   before the first and after the last. When `x` are the midpoints of
#'   histogram bins, as from [hist.weighted()], these are the bin edges.
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#' @examples
#' # Periodic smoothing over a week; the result has one value per break
#' breaks <- seq(0, 7, by = 0.05)
#' h <- hist.weighted(c(1, 1, 2, 5, 5, 5), rep(1, 6), breaks = breaks)
#' smooth <- ker.smooth.fft(h$mids, h$density, bandwidth = 0.5)
#' plot(breaks, smooth, type = "l", xlab = "Day of week", ylab = "Smoothed density")
#' @export
ker.smooth.fft <- function(x, z, bandwidth) {
  nx <- length(z)

  zz <- z[c(1:nx, seq(1, nx, 1))]

  dx <- x[2] - x[1]

  ker.val.x <- stats::dnorm((1:nx) * dx - dx / 2, 0, sd = bandwidth)
  ker.val.x <- c(ker.val.x, rev(ker.val.x))

  smoothed <- stats::fft(stats::fft(zz) * stats::fft(ker.val.x), inverse = TRUE) / nx / 2 * dx

  # With the kernel weights at half-step distances, element `i` of the
  # circular convolution is the smoothed value at `x[1] - dx / 2 + i * dx`.
  # Element `2 * nx`, the same as element 0 by periodicity, is therefore the
  # value at `x[1] - dx / 2`, and elements `1:nx` follow at steps of `dx`.
  Re(smoothed[c(2 * nx, 1:nx)])
}

#' @rdname ker.smooth.fft
#' @export
ker.smooth.1D.fft <- ker.smooth.fft

#' 2D Gaussian kernel smoothing by FFT, with reflected boundaries
#'
#' Smooths a matrix `z` sampled on an equally-spaced `x`/`y` grid with a
#' separable Gaussian kernel of standard deviations
#' `x.bandwidth`/`y.bandwidth`, using an FFT.
#'
#' Note that despite the shared "fft" in the name, this extends `z` by
#' *reflection* in both dimensions, whereas the 1D [ker.smooth.fft()] extends
#' periodically — appropriately so, since this one smooths the spatial
#' triggering kernel (not periodic) and that one smooths the weekly and
#' monthly background components (periodic). See [ker.smooth.fft()].
#'
#' @param x,y Numeric vectors of equally-spaced grid points for each
#'   dimension.
#' @param z Numeric matrix of values on the `x` by `y` grid
#'   (`nrow(z) == length(x)`, `ncol(z) == length(y)`).
#' @param x.bandwidth,y.bandwidth Standard deviations of the Gaussian
#'   smoothing kernel along `x` and `y`.
#'
#' @return A numeric matrix of dimension `(length(x)+1) x (length(y)+1)`:
#'   the smoothed values at the points halfway between successive `x` (and
#'   `y`), and half a grid step beyond the first and last. When `x` and `y`
#'   are the midpoints of histogram bins, as from [hist.weighted.2D()], these
#'   are the bin edges.
#'
#' @references D'Angelo et al. (2024); see [spHawkes-package].
#' @examples
#' # Smooth a single spike at the origin; the result, one more row and column
#' # than the input, is centred on the spike
#' g <- seq(-1000, 1000, by = 100)
#' z <- matrix(0, length(g), length(g))
#' z[11, 11] <- 1
#' smooth <- ker.smooth.2D.fft(g, g, z, x.bandwidth = 200, y.bandwidth = 200)
#' dim(smooth)
#' # The largest values sit half a step either side of the spike
#' edges <- c(g - 50, max(g) + 50)
#' peak <- which(smooth > max(smooth) * (1 - 1e-9), arr.ind = TRUE)
#' edges[unique(peak[, 1])]
#' @export
ker.smooth.2D.fft <- function(x, y, z, x.bandwidth, y.bandwidth) {
  nx <- nrow(z)
  ny <- ncol(z)

  zz <- z[c(1:nx, seq(nx, 1, -1)), c(1:ny, seq(ny, 1, -1))]

  dx <- x[2] - x[1]
  dy <- y[2] - y[1]

  ker.val.x <- stats::dnorm((1:nx) * dx - dx / 2, 0, sd = x.bandwidth)
  ker.val.y <- stats::dnorm((1:ny) * dy - dy / 2, 0, sd = y.bandwidth)
  ker.val.x <- c(ker.val.x, rev(ker.val.x))
  ker.val.y <- c(ker.val.y, rev(ker.val.y))

  smoothed <- stats::fft(stats::fft(zz) * stats::fft(ker.val.x %o% ker.val.y), inverse = TRUE) /
    nx / ny / 4 * dx * dy

  # As in ker.smooth.fft(): element `i` along each dimension is the value at
  # `x[1] - dx / 2 + i * dx`, so the values from `x[1] - dx / 2` onwards are
  # elements `2 * nx` (the same as element 0) and `1:nx`.
  Re(smoothed[c(2 * nx, 1:nx), c(2 * ny, 1:ny)])
}

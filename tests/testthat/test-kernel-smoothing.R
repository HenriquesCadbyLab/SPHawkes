test_that("simpson integrates a known function accurately", {
  x <- seq(0, pi, length.out = 101) # odd length required
  val <- sin(x)
  approx_integral <- simpson(val, dx = x[2] - x[1])
  expect_equal(approx_integral, 2, tolerance = 1e-4) # integral of sin(x) on [0,pi] = 2
})

test_that("simpson.2D integrates a known function accurately", {
  x <- seq(0, 1, length.out = 21)
  y <- seq(0, 1, length.out = 21)
  val <- outer(x, y, function(a, b) a * b) # integral over [0,1]x[0,1] = 0.25
  approx_integral <- simpson.2D(val, dx = x[2] - x[1], dy = y[2] - y[1])
  expect_equal(approx_integral, 0.25, tolerance = 1e-4)
})

test_that("simpson() and simpson.2D() reject an even number of points", {
  # Simpson's weights 1, 4, 2, ..., 2, 4, 1 only fit an odd number of points.
  # With an even number the 1D weights used to be recycled, giving a wrong
  # integral with only a recycling warning, and the 2D product failed with
  # "non-conformable arguments".
  x <- seq(0, 3, length.out = 4)
  expect_error(simpson(x^2, dx = 1), "odd number of equally spaced points")
  expect_error(simpson(c(1, 2)), "`value` has 2")
  expect_error(simpson.2D(matrix(1, 4, 5)), "`value` is 4 x 5")
  expect_error(simpson.2D(matrix(1, 5, 4)), "`value` is 5 x 4")

  # Odd counts are unaffected (and exact for a quadratic).
  x <- seq(0, 3, length.out = 5)
  expect_equal(simpson(x^2, dx = x[2] - x[1]), 9)
  expect_equal(simpson.2D(outer(x, x), dx = 0.75, dy = 0.75), 4.5^2)
})

test_that("hist.weighted matches base hist() with unit weights", {
  set.seed(1)
  x <- runif(500, 0, 10)
  breaks <- seq(0, 10, 1)
  base_h <- hist(x, breaks = breaks, plot = FALSE)
  weighted_h <- hist.weighted(x, rep(1, length(x)), breaks = breaks)
  expect_equal(as.numeric(weighted_h$counts), base_h$counts, tolerance = 1e-8)
  expect_equal(weighted_h$mids, base_h$mids, tolerance = 1e-8)
})

test_that("hist.weighted.1D is an alias for hist.weighted", {
  expect_identical(hist.weighted.1D, hist.weighted)
})

test_that("hist.weighted() and hist.weighted.2D() stop on values outside the breaks", {
  # The original wrapped this check in `try()`: it printed the error and
  # carried on, and the out-of-range value then surfaced as a recycling
  # warning or a length mismatch further down.
  expect_error(
    hist.weighted(c(-1, 5), c(1, 1), breaks = seq(0, 10, 1)),
    "`breaks` (0 to 10) do not span the range of `x` (-1 to 5)", fixed = TRUE
  )
  expect_error(
    hist.weighted.2D(c(0.5, 3), c(0.5, 0.5), c(1, 1), x.breaks = 0:2, y.breaks = 0:2),
    "do not span the range of `x` (0.5 to 3)", fixed = TRUE
  )
  # Values exactly on the outer breaks are still binned.
  h <- hist.weighted(c(0, 10), c(1, 1), breaks = seq(0, 10, 1))
  expect_equal(sum(h$counts), 2)
})

test_that("hist.weighted.2D matches base table-based binning with unit weights", {
  set.seed(2)
  x <- runif(300, 0, 5)
  y <- runif(300, 0, 5)
  xb <- seq(0, 5, 1)
  yb <- seq(0, 5, 1)
  out <- hist.weighted.2D(x, y, rep(1, length(x)), x.breaks = xb, y.breaks = yb)

  ix <- cut(x, breaks = xb, include.lowest = TRUE)
  iy <- cut(y, breaks = yb, include.lowest = TRUE)
  expected_counts <- table(ix, iy)

  expect_equal(sum(out$counts), sum(expected_counts))
  expect_equal(dim(out$counts), dim(expected_counts))
})

test_that("ker.smooth.conv and ker.smooth.fft both smooth sensibly", {
  # The two are alternate implementations (direct convolution vs FFT) and
  # are not guaranteed to be bit-identical (stats::convolve() correlates
  # rather than convolves by default), but both should smooth a noisy bump
  # into something close to the underlying signal.
  x <- seq(0, 10, 0.1)
  set.seed(3)
  truth <- dnorm(x, 5, 1)
  y <- truth + rnorm(length(x), 0, 0.01)

  conv_result <- ker.smooth.conv(x, y, bandwidth = 0.5)
  fft_result <- ker.smooth.fft(x, y, bandwidth = 0.5)

  expect_length(conv_result, length(x) + 1)
  expect_length(fft_result, length(x) + 1)
  expect_gt(stats::cor(conv_result[1:length(x)], truth), 0.99)
  expect_gt(stats::cor(fft_result[1:length(x)], truth), 0.99)
})

test_that("ker.smooth.fft() and ker.smooth.2D.fft() return values at the bin edges", {
  # A single nonzero bin smooths to a peak centred on that bin: the two
  # largest values are at its edges, and equal.
  breaks <- seq(0, 7, by = 0.05)
  mids <- head(breaks, -1) + 0.025
  z <- numeric(length(mids))
  z[40] <- 1
  s <- ker.smooth.fft(mids, z, bandwidth = 0.3)
  expect_length(s, length(breaks))
  top <- order(s, decreasing = TRUE)[1:2]
  expect_equal(sort(breaks[top]), breaks[40:41])
  expect_equal(s[top[1]], s[top[2]])
  # Periodic: the first and last edges are the same point of the cycle.
  expect_equal(s[1], s[length(s)])

  # A spike at the origin smooths to a surface symmetric about it.
  g <- seq(-1000, 1000, by = 100)
  z2 <- matrix(0, length(g), length(g))
  z2[11, 11] <- 1
  s2 <- ker.smooth.2D.fft(g, g, z2, x.bandwidth = 200, y.bandwidth = 200)
  expect_equal(dim(s2), c(22, 22))
  expect_equal(s2, s2[nrow(s2):1, ])
  expect_equal(s2, s2[, ncol(s2):1])
})

test_that("ker.smooth.1D.fft is an alias for ker.smooth.fft", {
  expect_identical(ker.smooth.1D.fft, ker.smooth.fft)
})

test_that("ker.smooth.fft smoothing preserves the total mass approximately", {
  x <- seq(0, 10, 0.1)
  y <- rep(0, length(x))
  y[50] <- 1 # a single spike
  smoothed <- ker.smooth.fft(x, y, bandwidth = 0.5)
  # Convolving with a Gaussian kernel (integrating to 1) should preserve
  # the total mass of the original series: sum(y) * dx.
  dx <- x[2] - x[1]
  expect_equal(sum(smoothed) * dx, sum(y) * dx, tolerance = 0.05 * dx)
})

test_that("ker.smooth.2D.fft returns a matrix of the expected dimensions", {
  x <- seq(0, 10, 1)
  y <- seq(0, 10, 1)
  z <- matrix(0, nrow = length(x), ncol = length(y))
  z[6, 6] <- 1
  smoothed <- ker.smooth.2D.fft(x, y, z, x.bandwidth = 1, y.bandwidth = 1)
  expect_equal(dim(smoothed), c(length(x) + 1, length(y) + 1))
})

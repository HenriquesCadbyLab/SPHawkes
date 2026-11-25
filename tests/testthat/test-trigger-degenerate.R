# Regression tests for triggering-kernel re-estimation when most offsets are
# unobservable.
#
# The spatial kernel is estimated by dividing an observed density by how often
# each offset could have been observed at all. For a study region smaller than
# the kernel's extent, that denominator is zero over most of the grid. The
# numerator is zero there too, and `0 / 0` gave `NaN`, which the FFT smoother
# then spread across the whole kernel: an all-`NA` kernel, `A = Inf`, and a fit
# that died with "missing value where TRUE/FALSE needed" from the convergence
# check. The synthetic region is 1000 x 1000; the tests that need that regime
# set a 4000 kernel extent explicitly (the default is 1000).

test_that(".ratio_or_zero treats an unusable denominator as no information", {
  expect_equal(.ratio_or_zero(c(1, 2, 0), c(2, 4, 0)), c(0.5, 0.5, 0))
  expect_equal(.ratio_or_zero(c(1, 1, 1), c(0, -1, NA)), c(0, 0, 0))
  expect_equal(.ratio_or_zero(c(1, 1), c(Inf, 2)), c(0, 0.5))

  # A non-finite numerator is dropped rather than propagated.
  expect_equal(.ratio_or_zero(c(NaN, Inf, 3), c(1, 1, 1)), c(0, 0, 3))

  # Shape is preserved, which matters because the spatial kernel is a matrix.
  m <- matrix(c(1, 2, 3, 4), 2)
  d <- matrix(c(1, 0, 2, 0), 2)
  r <- .ratio_or_zero(m, d)
  expect_equal(dim(r), c(2L, 2L))
  expect_equal(r, matrix(c(1, 0, 1.5, 0), 2))

  expect_equal(.ratio_or_zero(numeric(0), numeric(0)), numeric(0))
})

test_that(".normalise_kernel fails clearly instead of returning NaN", {
  expect_error(.normalise_kernel(c(1, 2), 0, "spatial"), "cannot be normalised")
  expect_error(.normalise_kernel(c(1, 2), NaN, "temporal"), "cannot be normalised")
  expect_error(.normalise_kernel(c(1, 2), -1, "spatial"), "cannot be normalised")
  expect_equal(.normalise_kernel(c(1, 2), 2, "spatial"), c(0.5, 1))
})

test_that("a fit on a region smaller than the kernel extent still produces finite estimates", {
  synth <- make_synthetic_data(n = 60, seed = 7)

  w <- expect_warning(
    fit <- fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 4000),
      n_iter = 1, verbose = FALSE, seed = 1
    ),
    "study region spans only"
  )
  # The suggested extent must be one trigger_spatial_init() accepts with the
  # default step: on this 1000-unit square it used to suggest 500, which
  # gives a 6-point grid (500 / 200 = 2.5 steps).
  suggested <- as.numeric(sub(".*x_range = ([0-9.]+).*", "\\1", conditionMessage(w)))
  expect_equal(suggested, 400)
  expect_no_error(trigger_spatial_init(x_range = suggested))

  expect_s3_class(fit, "semihawkes_fit")
  expect_true(is.finite(fit$mu) && fit$mu > 0)
  expect_true(is.finite(fit$A) && fit$A > 0)
  expect_true(is.finite(fit$llik))
  expect_false(anyNA(fit$bgprobs))
  # The kernels themselves must be usable, not silently all-NA.
  expect_false(anyNA(fit$spatial_trigger$basevalue))
  expect_false(anyNA(fit$temporal_trigger$basevalue))
})

test_that("trigger_spatial_args reaches trigger_spatial_init, and a matched extent fits cleanly", {
  synth <- make_synthetic_data(n = 60, seed = 7)

  expect_silent(
    fit <- fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 500, x_step = 50),
      trigger_window_dist = 500,
      n_iter = 1, verbose = FALSE, seed = 1
    )
  )

  expect_equal(max(fit$spatial_trigger$x), 500)
  expect_false(anyNA(fit$spatial_trigger$basevalue))
  expect_true(is.finite(fit$A) && fit$A > 0)
})

test_that("a pair window reaching past the kernel grid is rejected", {
  synth <- make_synthetic_data(n = 60, seed = 7)

  expect_error(
    fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 500, x_step = 50),
      trigger_window_dist = 1000,
      n_iter = 1, verbose = FALSE, seed = 1
    ),
    "larger than the spatial triggering kernel"
  )
})

test_that("an unset pair window follows the kernel grid's reach", {
  synth <- make_synthetic_data(n = 60, seed = 7)

  expect_silent(
    fit <- fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 500, x_step = 50),
      n_iter = 1, verbose = FALSE, seed = 1
    )
  )
  expect_equal(max(fit$spatial_trigger$x), 500)
  expect_true(is.finite(fit$A) && fit$A > 0)
})

test_that("a pair window longer than the temporal kernel's cutoff is rejected", {
  # With trigger_window_days = 120 and the 90-day kernel, hist.weighted()
  # used to print "'breaks' do not span range of 'x'" from inside `try()`
  # and the fit then failed with "length mismatch in convolution".
  synth <- make_synthetic_data(n = 60, seed = 7)

  expect_error(
    fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 500, x_step = 50),
      trigger_window_dist = 500, trigger_window_days = 120,
      n_iter = 1, verbose = FALSE, seed = 1
    ),
    "`trigger_window_days` (120) is longer than the temporal triggering kernel's cutoff (90 days)",
    fixed = TRUE
  )
  expect_error(
    fit_semihawkes(
      synth$points, synth$region, date_col = "event_date",
      cell_size = 250, cache_dir = tempfile("sphawkes-degen-"),
      bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
      trigger_spatial_args = list(x_range = 500, x_step = 50),
      trigger_window_dist = 500, trigger_window_days = "90",
      n_iter = 1, verbose = FALSE, seed = 1
    ),
    "must be a single number"
  )
})

test_that("bg_integral returns values in (0, 1]", {
  d <- make_synthetic_data(n = 100)
  ev <- prepare_events(d$points, d$region, "event_date")
  ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 50)

  vals <- bg_integral(ev$events$coorx, ev$events$coory, ev$events$bandwidth,
                       ev$region, n_mc = 500, seed = 42)
  expect_length(vals, nrow(ev$events))
  expect_true(all(vals > 0 & vals <= 1))
})

test_that("bg_integral gives values near 1 for points well inside a moderate region", {
  # Region-to-kernel area ratio comparable to the original city/bandwidth
  # scale (~1000x), so the Monte Carlo estimator has enough samples near
  # the event to be stable.
  square <- sf::st_polygon(list(matrix(
    c(-1500, -1500, 1500, -1500, 1500, 1500, -1500, 1500, -1500, -1500),
    ncol = 2, byrow = TRUE
  )))
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = 32717))
  val <- bg_integral(0, 0, 100, region, n_mc = 5000, seed = 1)
  expect_gt(val, 0.85)
})

test_that("add_bg_integral requires bandwidth to already be present", {
  d <- make_synthetic_data(n = 50)
  ev <- prepare_events(d$points, d$region, "event_date")
  expect_error(add_bg_integral(ev), "bandwidth")
})

test_that("add_bg_integral adds bg.integral and keeps row count consistent", {
  d <- make_synthetic_data(n = 100)
  ev <- prepare_events(d$points, d$region, "event_date")
  ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 50)
  ev <- add_bg_integral(ev, n_mc = 500, seed = 7)

  expect_true("bg.integral" %in% names(ev$events))
  expect_equal(ev$n_kept, nrow(ev$events))
  expect_true(all(ev$events$bg.integral > 0 & ev$events$bg.integral <= 1))
})

test_that("a seeded bg_integral() is reproducible and leaves the caller's RNG state alone", {
  square <- sf::st_polygon(list(matrix(
    c(0, 0, 1000, 0, 1000, 1000, 0, 1000, 0, 0),
    ncol = 2, byrow = TRUE
  )))
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = 32717))

  a <- bg_integral(c(100, 500), c(100, 500), c(200, 200), region,
                   n_mc = 500, seed = 3)
  b <- bg_integral(c(100, 500), c(100, 500), c(200, 200), region,
                   n_mc = 500, seed = 3)
  expect_identical(a, b)

  set.seed(42)
  expected_next <- stats::runif(3)
  set.seed(42)
  bg_integral(100, 100, 200, region, n_mc = 500, seed = 3)
  expect_identical(stats::runif(3), expected_next)
})

test_that("bg_integral(restore_rng = FALSE) leaves the generator in the seeded state", {
  square <- sf::st_polygon(list(matrix(
    c(0, 0, 1000, 0, 1000, 1000, 0, 1000, 0, 0),
    ncol = 2, byrow = TRUE
  )))
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = 32717))

  set.seed(42)
  restored_next <- stats::runif(3)

  # Whatever the state beforehand, the state afterwards is fixed by `seed`
  set.seed(42)
  bg_integral(100, 100, 200, region, n_mc = 500, seed = 3, restore_rng = FALSE)
  after_a <- stats::runif(3)
  set.seed(7)
  bg_integral(100, 100, 200, region, n_mc = 500, seed = 3, restore_rng = FALSE)
  after_b <- stats::runif(3)

  expect_identical(after_a, after_b)
  expect_false(identical(after_a, restored_next))

  # The integral itself does not depend on restore_rng
  expect_identical(
    bg_integral(100, 100, 200, region, n_mc = 500, seed = 3, restore_rng = FALSE),
    bg_integral(100, 100, 200, region, n_mc = 500, seed = 3)
  )
  expect_error(
    bg_integral(100, 100, 200, region, n_mc = 500, seed = 3, restore_rng = NA),
    "restore_rng"
  )
})

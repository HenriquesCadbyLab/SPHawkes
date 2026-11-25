test_that("adaptive_bandwidth returns min_bandwidth for isolated points", {
  set.seed(1)
  x <- c(0, 1000, 2000, 3000, 4000, 5000, 6000) # far apart
  y <- rep(0, length(x))
  bw <- adaptive_bandwidth(x, y, k = 2, min_dist = 1, min_bandwidth = 250)
  expect_true(all(bw >= 250))
})

test_that("adaptive_bandwidth ignores near-duplicate points via min_dist", {
  # 6 points at the origin (duplicates), 1 point far away.
  x <- c(rep(0, 6), 500)
  y <- c(rep(0, 6), 0)
  bw <- adaptive_bandwidth(x, y, k = 1, min_dist = 150, min_bandwidth = 10)
  # The lone far point's nearest neighbour beyond min_dist should be ~500.
  expect_equal(bw[7], 500, tolerance = 1e-6)
})

test_that("adaptive_bandwidth warns and falls back when too few neighbours exist", {
  x <- c(0, 10, 4000)
  y <- c(0, 0, 0)
  expect_warning(
    bw <- adaptive_bandwidth(x, y, k = 2, min_dist = 150, min_bandwidth = 99),
    "neighbours"
  )
  expect_true(all(bw >= 99))
})

test_that("add_adaptive_bandwidth requires a semihawkes_events object", {
  expect_error(add_adaptive_bandwidth(list(events = data.frame())),
               "semihawkes_events")
})

test_that("add_adaptive_bandwidth adds a bandwidth column", {
  d <- make_synthetic_data(n = 150)
  ev <- prepare_events(d$points, d$region, "event_date")
  ev2 <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 250)
  expect_true("bandwidth" %in% names(ev2$events))
  expect_true(all(ev2$events$bandwidth >= 250))
  expect_equal(nrow(ev2$events), nrow(ev$events))
})

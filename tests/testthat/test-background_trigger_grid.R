make_grid_fit_inputs <- function(n = 250, seed = 21) {
  synth <- make_synthetic_data(n = n, seed = seed)

  cell_size <- 150
  city <- sf::st_buffer(synth$region, dist = 4 * cell_size)

  fit <- fit_semihawkes(
    points = synth$points,
    region = city,
    date_col = "event_date",
    cell_size = cell_size,
    cache_dir = tempfile("spHawkes-grid-marks-"),
    bandwidth_args = list(
      k = 3,
      min_dist = 10,
      min_bandwidth = 50
    ),
    n_iter = 0,
    verbose = FALSE
  )

  list(fit = fit, city = city)
}

test_that("background_trigger_grid returns correctly shaped output", {
  inputs <- make_grid_fit_inputs()

  out <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150,
                                 days = 0:20)

  expect_type(out, "list")
  expect_named(out, c("days", "x", "y", "inside", "snapshots"))
  expect_equal(out$days, 0:20)
  expect_length(out$snapshots, length(0:20))
  expect_length(out$inside, length(out$x) * length(out$y))

  # Every daily snapshot must be a grid-shaped matrix for both components
  first <- out$snapshots[[1]]
  expect_named(first, c("bg", "tr"))
  expect_equal(dim(first$bg), c(length(out$x), length(out$y)))
  expect_equal(dim(first$tr), c(length(out$x), length(out$y)))
})

test_that("background_trigger_grid's intensities are finite and non-negative", {
  inputs <- make_grid_fit_inputs()

  out <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150,
                                 days = 0:20)

  all_bg <- unlist(lapply(out$snapshots, function(s) s$bg))
  all_tr <- unlist(lapply(out$snapshots, function(s) s$tr))

  expect_true(all(is.finite(all_bg)))
  expect_true(all(is.finite(all_tr)))
  expect_true(all(all_bg >= 0))
  expect_true(all(all_tr >= 0))
})

test_that("background_trigger_grid sets cells outside the study region to zero", {
  inputs <- make_grid_fit_inputs()

  out <- background_trigger_grid(
    inputs$fit,
    inputs$city,
    cell_size = 150,
    days = 0:5
  )

  outside <- which(!out$inside)

  # The padded synthetic city guarantees that the grid extends beyond the
  # study region, so there must be cells outside it.
  expect_true(length(outside) > 0)

  for (snap in out$snapshots) {
    expect_true(all(snap$bg[outside] == 0))
    expect_true(all(snap$tr[outside] == 0))
  }
})

test_that("background_trigger_grid's triggering component reacts to t.lag", {
  inputs <- make_grid_fit_inputs()

  # With t.lag = 0, no past event should ever be close enough in time to
  # contribute -- so the triggering component must be identically zero at
  # every grid cell, on every day.
  out_no_lag <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150,
                                        days = 0:10, t.lag = 0)

  all_tr_no_lag <- unlist(lapply(out_no_lag$snapshots, function(s) s$tr))
  expect_true(all(all_tr_no_lag == 0))
})

test_that("background_trigger_grid's background component matches a manual calculation at one cell", {
  inputs <- make_grid_fit_inputs()
  fit <- inputs$fit

  # Days in different weeks and months, so a weekly or monthly term evaluated
  # at the wrong point in the calendar cannot match by coincidence.
  days <- c(5, 40, 100, 200)
  out <- background_trigger_grid(fit, inputs$city, cell_size = 150, days = days)

  inside_idx <- which(out$inside)
  skip_if(length(inside_idx) == 0, "no inside cells to check")
  check_idx <- inside_idx[1]

  grid <- expand.grid(x = out$x, y = out$y)
  dates <- fit$events$start_date + days
  wday <- as.POSIXlt(dates)$wday
  month <- as.POSIXlt(dates)$mon

  expected_bg <- fit$mu * fit$trend$fun(days) * fit$weekly$fun(wday) *
    fit$monthly$fun(month) *
    fit$background$fun(grid$x[check_idx], grid$y[check_idx])

  actual_bg <- vapply(out$snapshots, function(s) s$bg[check_idx], numeric(1))
  # Intensities here are of order 1e-13, below the point where
  # `expect_equal()` switches from relative to absolute tolerance, so compare
  # the ratio to 1 rather than the raw values.
  expect_true(all(expected_bg > 0))
  expect_equal(actual_bg / expected_bg, rep(1, length(days)), tolerance = 1e-10)
})

test_that("background_trigger_grid evaluates weekly/monthly terms by calendar date", {
  inputs <- make_grid_fit_inputs()
  fit <- inputs$fit

  # This check is only sensitive if the study does not start on a Sunday
  # (day offset 0 would then coincide with weekday 0 by luck) and the fitted
  # weekly term is not flat.
  expect_false(as.POSIXlt(fit$events$start_date)$wday == 0)
  expect_gt(diff(range(fit$weekly$basevalue)), 1e-6)

  days <- 0:13
  out <- background_trigger_grid(fit, inputs$city, cell_size = 150, days = days)
  totals <- vapply(out$snapshots, function(s) sum(s$bg), numeric(1))

  dates <- fit$events$start_date + days
  calendar_factor <- fit$trend$fun(days) *
    fit$weekly$fun(as.POSIXlt(dates)$wday) *
    fit$monthly$fun(as.POSIXlt(dates)$mon)
  raw_day_factor <- fit$trend$fun(days) * fit$weekly$fun(days) *
    fit$monthly$fun(days)

  # Grid totals share the spatial surface across days, so they must be
  # proportional to the calendar-based temporal factor ... (compared after
  # rescaling, since the raw totals are too small for a relative tolerance)
  ratio <- totals / calendar_factor
  expect_equal(ratio / ratio[1], rep(1, length(ratio)), tolerance = 1e-10)
  # ... and not to the raw-day-offset one.
  raw_ratio <- totals / raw_day_factor
  expect_gt(diff(range(raw_ratio)) / mean(raw_ratio), 1e-3)
})

test_that("background_trigger_grid defaults `days` to 0:TT from the fit's events", {
  inputs <- make_grid_fit_inputs()

  out <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150)

  expected_TT <- max(inputs$fit$events$events$days) + 1
  expect_equal(out$days, 0:expected_TT)
})

test_that("background_trigger_grid(nprocs = 2) does not leave its stopped cluster registered", {
  # It used to stop its cluster on exit but leave it registered as the
  # foreach backend, so the next `%dopar%` in the session failed with
  # "invalid connection".
  skip_on_cran()
  skip_if_not_installed("doParallel")
  inputs <- make_grid_fit_inputs()

  parallel_out <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150,
                                          days = 0:3, nprocs = 2)
  expect_equal(foreach::getDoParName(), "doSEQ")
  `%dopar%` <- foreach::`%dopar%`
  expect_equal(foreach::foreach(i = 1:3, .combine = c) %dopar% i^2, c(1, 4, 9))

  serial_out <- background_trigger_grid(inputs$fit, inputs$city, cell_size = 150,
                                        days = 0:3, nprocs = 1)
  expect_equal(parallel_out$snapshots, serial_out$snapshots)
})

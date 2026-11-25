make_ready_events <- function(n = 60, seed = 7) {
  synth <- make_synthetic_data(n = n, seed = seed)
  ev <- prepare_events(synth$points, synth$region, date_col = "event_date",
                        start_date = "2020-01-01")
  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
  ev
}

test_that("background_spatial requires cell_size explicitly", {
  ev <- make_ready_events()
  expect_error(background_spatial(ev), "cell_size")
  expect_error(background_spatial(ev, cell_size = -10), "cell_size")
})

test_that("background_spatial requires bandwidth/bg.integral to already exist", {
  synth <- make_synthetic_data(n = 30, seed = 3)
  ev <- prepare_events(synth$points, synth$region, date_col = "event_date",
                        start_date = "2020-01-01")
  expect_error(background_spatial(ev, cell_size = 100), "bandwidth")

  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  expect_error(background_spatial(ev, cell_size = 100), "bg.integral")
})

test_that("background_spatial returns a mean-1 (over in-region cells) grid", {
  ev <- make_ready_events()
  sb <- background_spatial(ev, cell_size = 100)

  expect_s3_class(sb, "semihawkes_spatial_base")
  expect_equal(dim(sb$basevalue), c(length(sb$x), length(sb$y)))
  expect_equal(dim(sb$marks), dim(sb$basevalue))
  expect_equal(mean(sb$basevalue[sb$marks > 0]), 1, tolerance = 1e-8)
  expect_true(all(sb$marks %in% c(0L, 1L)))
})

test_that("background_spatial chunking doesn't change the result", {
  ev <- make_ready_events(n = 45)
  sb_big_chunks <- background_spatial(ev, cell_size = 120, chunk_size = 1000)
  sb_small_chunks <- background_spatial(ev, cell_size = 120, chunk_size = 5)
  expect_equal(sb_big_chunks$basevalue, sb_small_chunks$basevalue, tolerance = 1e-10)
})

test_that("background_spatial matches a brute-force per-event loop", {
  ev <- make_ready_events(n = 25)
  resp_weight <- stats::runif(nrow(ev$events), 0.3, 1)
  sb <- background_spatial(ev, cell_size = 150, resp_weight = resp_weight,
                            chunk_size = 4)

  basex <- sb$x %o% rep(1, length(sb$y))
  basey <- rep(1, length(sb$x)) %o% sb$y
  brute <- matrix(0, nrow = length(sb$x), ncol = length(sb$y))
  for (i in seq_len(nrow(ev$events))) {
    bgsmoother <- stats::dnorm(basex, ev$events$coorx[i], ev$events$bandwidth[i]) *
      stats::dnorm(basey, ev$events$coory[i], ev$events$bandwidth[i]) /
      ev$events$bg.integral[i]
    brute <- brute + resp_weight[i] * bgsmoother
  }
  brute <- brute / mean(brute[sb$marks > 0])

  expect_equal(sb$basevalue, brute, tolerance = 1e-8)
})

test_that("background_spatial's fun bilinearly interpolates and is NA outside the grid", {
  ev <- make_ready_events()
  sb <- background_spatial(ev, cell_size = 100)

  # a grid vertex should be reproduced exactly
  val <- sb$fun(sb$x[3], sb$y[4])
  expect_equal(as.numeric(val), sb$basevalue[3, 4], tolerance = 1e-10)

  # far outside the grid -> NA
  outside <- sb$fun(min(sb$x) - 1e6, min(sb$y) - 1e6)
  expect_true(is.na(outside))
})

test_that("background_intensity combines temporal + spatial base values", {
  ev <- make_ready_events(n = 80)
  wb <- temporal_base_weekly(ev)
  mb <- temporal_base_monthly(ev)
  tb <- temporal_base_trend(ev)
  sb <- background_spatial(ev, cell_size = 120)

  bi <- background_intensity(ev, trend = tb, weekly = wb, monthly = mb, spatial = sb)

  expect_length(bi$mub_events, nrow(ev$events))
  expect_length(bi$bg_at_events, nrow(ev$events))
  expect_equal(dim(bi$mub_grid), dim(sb$basevalue))
  expect_true(is.numeric(bi$bg_at_all) && length(bi$bg_at_all) == 1)
  expect_true(bi$bg_at_all > 0)

  # bg_at_events should equal the product of the four factors, event by event
  dates <- ev$start_date + ev$events$days
  dt <- as.POSIXlt(dates)
  expected <- tb$fun(ev$events$days) * wb$fun(dt$wday) *
    mb$fun(dt$mon) * bi$mub_events
  expect_equal(bi$bg_at_events, expected, tolerance = 1e-10)
})

test_that("background_intensity validates its inputs' classes", {
  ev <- make_ready_events()
  wb <- temporal_base_weekly(ev)
  mb <- temporal_base_monthly(ev)
  tb <- temporal_base_trend(ev)
  sb <- background_spatial(ev, cell_size = 120)

  expect_error(background_intensity(ev, trend = "nope", weekly = wb, monthly = mb,
                                     spatial = sb),
               "semihawkes_temporal_base")
  expect_error(background_intensity(ev, trend = tb, weekly = wb, monthly = mb,
                                     spatial = "nope"),
               "semihawkes_spatial_base")
})

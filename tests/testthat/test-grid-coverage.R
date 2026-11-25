# Regression tests for the background grid covering the study area.
#
# `seq(from, to, by)` stops at or before `to`, so a grid built with it fell
# short of the bounding box whenever `cell_size` did not divide the bbox
# exactly. Events in the uncovered strip then had no background value, and the
# fit died inside an interpolation call with a message that named neither
# `cell_size` nor the grid. The synthetic region is 1000 x 1000, so the cell
# sizes used below deliberately do not divide it.

ready_events <- function(n = 60, seed = 7) {
  synth <- make_synthetic_data(n = n, seed = seed)
  ev <- prepare_events(synth$points, synth$region, date_col = "event_date",
                       start_date = "2020-01-01")
  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  add_bg_integral(ev, n_mc = 500, seed = 1)
}

test_that(".cover_seq always reaches its endpoint", {
  # Where the step divides the span exactly, the grid is unchanged.
  expect_equal(.cover_seq(0, 1000, 250), seq(0, 1000, by = 250))
  expect_equal(.cover_seq(0, 8000, 500), seq(0, 8000, by = 500))

  for (by in c(97, 150, 220, 7, pi)) {
    g <- .cover_seq(0, 1000, by)
    expect_gte(max(g), 1000)                       # covers the endpoint
    expect_lt(max(g) - 1000, by)                   # by less than one cell
    expect_equal(diff(g), rep(by, length(g) - 1))  # still regular
    expect_equal(g[1], 0)
  }

  # A step larger than the span still gives a usable, covering grid.
  g <- .cover_seq(0, 10, 1000)
  expect_gte(length(g), 2)
  expect_gte(max(g), 10)
})

test_that("the background grid covers the bounding box for any cell_size", {
  ev <- ready_events()
  bbox <- sf::st_bbox(ev$region)

  for (cs in c(97, 150, 220)) {
    sb <- background_spatial(ev, cell_size = cs)
    expect_gte(max(sb$x), bbox[["xmax"]])
    expect_gte(max(sb$y), bbox[["ymax"]])
    # No event is left without a background value.
    expect_false(anyNA(sb$fun(ev$events$coorx, ev$events$coory)))
  }
})

test_that("the background integral does not depend on how cell_size divides the region", {
  ev <- ready_events()
  tb <- temporal_base_trend(ev)
  wb <- temporal_base_weekly(ev)
  mb <- temporal_base_monthly(ev)

  bg_at_all <- vapply(c(97, 125, 150, 220, 250), function(cs) {
    background_intensity(ev, tb, wb, mb,
                         background_spatial(ev, cell_size = cs))$bg_at_all
  }, numeric(1))

  # The background shape is normalised to mean 1 over the in-region cells, so
  # its spatial integral is the region's area at any grid resolution.
  expect_equal(max(bg_at_all) / min(bg_at_all), 1, tolerance = 1e-8)
})

test_that(".grid_cell_area tiles the region exactly", {
  ev <- ready_events()
  for (cs in c(97, 150, 220)) {
    sb <- background_spatial(ev, cell_size = cs)
    expect_equal(sum(sb$marks > 0) * .grid_cell_area(sb, ev$region),
                 as.numeric(sf::st_area(ev$region)))
  }
})

test_that("background_intensity rejects a grid that does not cover the events", {
  ev <- ready_events()

  # Build a background surface on a sub-region, so its grid genuinely does not
  # reach the events in the rest of the study area.
  synth <- make_synthetic_data(n = 60, seed = 7)
  small <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(matrix(
    c(0, 0, 500, 0, 500, 500, 0, 500, 0, 0), ncol = 2, byrow = TRUE
  ))), crs = 32717))
  ev_small <- prepare_events(synth$points, small, date_col = "event_date",
                             start_date = "2020-01-01")
  ev_small <- add_adaptive_bandwidth(ev_small, k = 3, min_dist = 10,
                                     min_bandwidth = 50)
  ev_small <- add_bg_integral(ev_small, n_mc = 500, seed = 1)
  sb_small <- background_spatial(ev_small, cell_size = 100)

  expect_error(
    background_intensity(ev, temporal_base_trend(ev), temporal_base_weekly(ev),
                         temporal_base_monthly(ev), sb_small),
    "outside the spatial background grid"
  )
})

# The two end-to-end tests below use a region much larger than the default
# triggering-kernel extent (+/-1000 units from `trigger_spatial_init()`). On a
# much smaller region the edge-correction denominator is zero over most of that
# grid, which is a separate problem from the one these tests cover.
wide_data <- function(n = 80, seed = 3) {
  set.seed(seed)
  square <- sf::st_polygon(list(matrix(
    c(0, 0, 9000, 0, 9000, 9000, 0, 9000, 0, 0), ncol = 2, byrow = TRUE
  )))
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = 32717))
  pts <- sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 9000), y = stats::runif(n, 0, 9000),
               event_date = as.Date("2020-01-01") + sample(0:364, n, replace = TRUE)),
    coords = c("x", "y"), crs = 32717
  )
  list(points = pts, region = region)
}

test_that("fit_semihawkes runs with a cell_size that does not divide the region", {
  d <- wide_data()
  # 9000 / 1300 is not a whole number, so the old grid stopped short of the
  # bounding box and events in the last strip had no background value.
  fit <- fit_semihawkes(
    d$points, d$region, date_col = "event_date",
    cell_size = 1300, cache_dir = tempfile("sphawkes-grid-test-"),
    bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 300),
    n_iter = 1, verbose = FALSE, seed = 1
  )
  expect_s3_class(fit, "semihawkes_fit")
  expect_false(anyNA(fit$bgprobs))
  expect_true(is.finite(fit$mu) && fit$mu > 0)
  expect_true(is.finite(fit$A) && fit$A > 0)
  expect_true(is.finite(fit$llik))
})

test_that("fit_semihawkes is reproducible when given a seed", {
  d <- wide_data()
  run <- function() fit_semihawkes(
    d$points, d$region, date_col = "event_date",
    cell_size = 1500, cache_dir = tempfile("sphawkes-seed-test-"),
    bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 300),
    n_iter = 1, verbose = FALSE, seed = 42
  )
  a <- run(); b <- run()
  expect_equal(a$mu, b$mu)
  expect_equal(a$A, b$A)
  expect_equal(a$llik, b$llik)
})

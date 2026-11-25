fit <- make_gof_fit()

test_that("gof_spatial_residuals() validates its inputs", {
  expect_error(
    gof_spatial_residuals(fit, start_date = "2025-01-01"),
    "must define a non-empty window"
  )

  expect_error(
    gof_spatial_residuals(
      fit,
      start_date = "2024-01-01",
      end_date = "2023-01-01"
    ),
    "must define a non-empty window"
  )
})

test_that("gof_spatial_residuals() returns a coherent object", {
  gof <- gof_spatial_residuals(fit)

  expect_s3_class(gof, "semihawkes_gof_spatial")

  expect_named(
    gof$cells,
    c(
      "x",
      "y",
      "expected",
      "observed",
      "residual",
      "relative_residual",
      "standardized_residual"
    )
  )

  expect_true(all(gof$cells$expected >= 0))
  expect_true(all(gof$cells$observed >= 0))
  expect_true(all(gof$cells$observed == round(gof$cells$observed)))
  expect_true(all(gof$cells$residual ==
                    gof$cells$observed - gof$cells$expected))

  expect_equal(
    sum(gof$cells$observed),
    nrow(fit$events$events),
    tolerance = 1e-12
  )

  expect_equal(
    sum(gof$cells$expected),
    fit$lambda_at_all,
    tolerance = 0.05
  )
})

test_that("plot method returns a ggplot object", {
  gof <- gof_spatial_residuals(fit)

  p <- plot(gof)

  expect_s3_class(p, "ggplot")
})

test_that("plot method accepts custom limits", {
  gof <- gof_spatial_residuals(fit)

  p <- plot(gof, limits = c(-10, 10))

  expect_s3_class(p, "ggplot")
})

test_that("spatial residuals include relative and standardized versions", {
  gof <- gof_spatial_residuals(fit)

  expect_named(
    gof$cells,
    c(
      "x", "y", "expected", "observed", "residual",
      "relative_residual", "standardized_residual"
    )
  )

  expect_equal(
    gof$cells$residual,
    gof$cells$observed - gof$cells$expected,
    tolerance = 1e-12
  )

  expect_equal(
    gof$cells$relative_residual,
    (gof$cells$observed - gof$cells$expected) /
      gof$cells$expected,
    tolerance = 1e-12
  )

  expect_equal(
    gof$cells$standardized_residual,
    (gof$cells$observed - gof$cells$expected) /
      sqrt(gof$cells$expected),
    tolerance = 1e-12
  )
})

test_that("plot method supports residual types", {
  gof <- gof_spatial_residuals(fit)

  expect_s3_class(plot(gof, type = "raw"), "ggplot")
  expect_s3_class(plot(gof, type = "relative"), "ggplot")
  expect_s3_class(plot(gof, type = "standardized"), "ggplot")

  expect_error(
    plot(gof, type = "invalid"),
    "'arg' should be one of"
  )
})

test_that("the window includes its last day on both the observed and the expected side", {
  # Events on `end_date` used to be counted while the expected count stopped
  # at the start of that day, and a one-day window was an error.
  d <- fit$events$start_date
  one_day <- gof_spatial_residuals(fit, start_date = d, end_date = d)
  expect_equal(one_day$window_length, 1)
  expect_equal(one_day$end_date, d)
  expect_equal(sum(one_day$cells$observed), sum(fit$events$events$date == d))

  # Consecutive windows add up to the window covering both.
  a <- gof_spatial_residuals(fit, start_date = d, end_date = d + 9)
  b <- gof_spatial_residuals(fit, start_date = d + 10, end_date = d + 19)
  ab <- gof_spatial_residuals(fit, start_date = d, end_date = d + 19)
  expect_equal(ab$window_length, 20)
  expect_equal(sum(a$cells$observed) + sum(b$cells$observed), sum(ab$cells$observed))
  expect_equal(sum(a$cells$expected) + sum(b$cells$expected), sum(ab$cells$expected),
               tolerance = 1e-3)

  # The default window is the whole observation window, as before.
  whole <- gof_spatial_residuals(fit)
  expect_equal(c(whole$start_day, whole$end_day), c(0, fit$events$TT))
})

test_that("the residual plots draw cells at the grid's cell size", {
  # The plots took the cell size from the gap between the first two x (and y)
  # values of the cells in the region. When the second column of the grid has
  # no cell in the region, as for a detached piece of the region at its edge,
  # every tile came out two cells wide.
  gof <- gof_spatial_residuals(fit)
  expect_equal(gof$cell_size, fit$background$cell_size)

  columns <- sort(unique(gof$cells$x))
  gappy <- gof
  gappy$cells <- gof$cells[gof$cells$x != columns[2], ]

  tile_widths <- function(p) {
    tiles <- sf::st_geometry(p$layers[[1]]$data)
    vapply(tiles, function(g) {
      bb <- sf::st_bbox(g)
      bb[["xmax"]] - bb[["xmin"]]
    }, numeric(1))
  }
  expect_equal(max(tile_widths(plot(gappy))), fit$background$cell_size)

  gappy_periods <- gof_period_residuals(fit, period = "year")
  gappy_periods[[1]]$cells <- gappy_periods[[1]]$cells[gappy_periods[[1]]$cells$x != columns[2], ]
  expect_length(gappy_periods, 1)
  expect_equal(max(tile_widths(plot(gappy_periods))), fit$background$cell_size)

  # Residual objects made before `cell_size` was stored use the smallest gap.
  gappy$cell_size <- NULL
  expect_equal(max(tile_widths(plot(gappy))), fit$background$cell_size)
})

test_that("the expected counts match a sum over every event and every cell", {
  # gof_spatial_residuals() skips events and cells the triggering kernels
  # cannot reach; this adds up every event's triggering over every cell, as
  # it used to, for a window well inside the observation window.
  d <- fit$events$start_date
  gof <- gof_spatial_residuals(fit, start_date = d + 120, end_date = d + 139)

  ev <- fit$events$events
  bg <- fit$background
  marks <- fit$events$time.marks
  marks <- marks[marks >= 120 & marks <= 140]
  grid_x <- bg$x %o% rep(1, length(bg$y))
  grid_y <- rep(1, length(bg$x)) %o% bg$y
  cell_area <- .grid_cell_area(bg, fit$events$region)
  calendar <- .calendar_marks(marks, fit$events$start_date)
  background <- fit$mu * mean(fit$trend$fun(marks) * fit$weekly$fun(calendar$wday) *
                                fit$monthly$fun(calendar$month)) * 20 *
    bg$basevalue * bg$marks * cell_area
  triggering <- matrix(0, length(bg$x), length(bg$y))
  for (i in which(ev$days < 140)) {
    triggering <- triggering +
      mean(fit$temporal_trigger$fun(marks - ev$days[i])) * 20 *
      (bg$marks > 0) * fit$spatial_trigger$fun(grid_x - ev$coorx[i], grid_y - ev$coory[i]) *
      cell_area
  }
  brute <- (background + fit$A * triggering)[bg$marks > 0]

  expect_equal(gof$cells$expected, brute, tolerance = 1e-12)
})

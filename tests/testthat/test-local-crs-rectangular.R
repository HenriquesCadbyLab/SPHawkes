# A small study region that is (a) rectangular, so the background grid has
# different numbers of rows and columns, and (b) in a local engineering CRS
# like the bundled `crimes` data, for which GDAL/PROJ cannot build any
# coordinate transformation. Both used to break the spatial diagnostics and
# map plots; `sim_region` is square and projected, so it never exposed them.
make_rect_local_fit <- function(n = 150, seed = 3) {
  set.seed(seed)
  rect <- sf::st_polygon(list(matrix(
    c(0, 0, 2000, 0, 2000, 1000, 0, 1000, 0, 0),
    ncol = 2, byrow = TRUE
  )))
  crs <- sf::st_crs(crimes_region)
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(rect, crs = crs))
  pts <- sf::st_as_sf(
    data.frame(
      x = stats::runif(n, 0, 2000),
      y = stats::runif(n, 0, 1000),
      date = as.Date("2020-01-01") + sample(0:364, n, replace = TRUE)
    ),
    coords = c("x", "y"), crs = crs
  )
  fit_semihawkes(
    pts, region, "date",
    cell_size = 250,
    cache_dir = tempfile("spHawkes-rect-"),
    bandwidth_args = list(k = 3, min_dist = 10, min_bandwidth = 50),
    trigger_spatial_args = list(x_range = 1000, x_step = 100),
    trigger_window_dist = 1000,
    n_iter = 0,
    verbose = FALSE
  )
}

test_that("gof_spatial_residuals() works on a rectangular background grid", {
  fit <- make_rect_local_fit()
  expect_false(length(fit$background$x) == length(fit$background$y))

  res <- gof_spatial_residuals(fit)
  expect_s3_class(res, "semihawkes_gof_spatial")
  expect_equal(nrow(res$cells), sum(fit$background$marks > 0))
  expect_true(all(is.finite(res$cells$expected)))
  expect_true(all(res$cells$expected >= 0))

  expect_s3_class(gof_period_residuals(fit, period = "month"),
                  "semihawkes_gof_period")
})

test_that("map plots work for data in a local engineering CRS", {
  fit <- make_rect_local_fit()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot(fit, which = "background"))
  expect_no_error(print(plot(gof_spatial_residuals(fit))))
  expect_no_error(print(plot(gof_period_residuals(fit, period = "month"))))
})

test_that(".plot_safe_crs() drops only engineering CRSs", {
  expect_true(is.na(sf::st_crs(.plot_safe_crs(crimes_region))))
  expect_identical(sf::st_crs(.plot_safe_crs(sim_region)),
                   sf::st_crs(sim_region))
})

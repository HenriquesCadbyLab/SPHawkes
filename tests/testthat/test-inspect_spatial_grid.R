test_that("inspect_spatial_grid returns the expected grid", {
  # Rectangle: 1,000 x 600 metres
  region <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = 0, ymin = 0, xmax = 1000, ymax = 600),
      crs = sf::st_crs(32717)
    )
  )

  result <- inspect_spatial_grid(
    region,
    cell_size = 200,
    plot = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result$grid, "sf")
  expect_equal(result$width_km, 1)
  expect_equal(result$height_km, 0.6)
  expect_equal(result$n_cells, 15L)
  expect_equal(nrow(result$grid), 15L)
  expect_equal(sf::st_crs(result$grid), sf::st_crs(region))
})


test_that("inspect_spatial_grid rejects invalid cell sizes", {
  region <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = 0, ymin = 0, xmax = 1000, ymax = 600),
      crs = sf::st_crs(32717)
    )
  )

  expect_error(
    inspect_spatial_grid(region, cell_size = 0, plot = FALSE),
    "single positive number"
  )
})

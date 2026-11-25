## tests/testthat/test-background-plot-utils.R

## --- Helpers to build minimal synthetic inputs -----------------------------

make_test_region <- function(crs = 27700) {
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)
    ))),
    crs = crs
  )
}

make_test_background <- function() {
  x <- seq(-2, 12, length.out = 8)
  y <- seq(-2, 12, length.out = 8)
  basevalue <- outer(x, y, function(xi, yi) xi + yi)
  marks <- matrix(1, nrow = length(x), ncol = length(y))
  marks[1, ] <- 0  # simulate a "no data" edge row, as background_spatial() would
  list(x = x, y = y, basevalue = basevalue, marks = marks)
}

## --- background_to_grid_sf() ------------------------------------------------

test_that("background_to_grid_sf() builds an sf grid with matching CRS", {
  region <- make_test_region()
  background <- make_test_background()

  grid_sf <- background_to_grid_sf(background, region)

  expect_s3_class(grid_sf, "sf")
  expect_equal(sf::st_crs(grid_sf), sf::st_crs(region))
  expect_true(all(c("x", "y", "basevalue") %in% names(grid_sf)))
  expect_equal(nrow(grid_sf), length(background$x) * length(background$y))
})

test_that("background_to_grid_sf() sets marks<=0 cells to NA", {
  region <- make_test_region()
  background <- make_test_background()

  grid_sf <- background_to_grid_sf(background, region)

  na_rows <- grid_sf[is.na(grid_sf$basevalue), ]
  expect_true(all(na_rows$x %in% background$x[1]))
})

## --- crop_background_spatial() ---------------------------------------------

test_that("crop_background_spatial() returns an sf object clipped to region", {
  region <- make_test_region()
  background <- make_test_background()

  grid_cropped <- crop_background_spatial(background, region, crop = TRUE)

  expect_s3_class(grid_cropped, "sf")
  expect_true(all(sf::st_geometry_type(grid_cropped) %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("crop_background_spatial() output falls within region's bounding box", {
  region <- make_test_region()
  background <- make_test_background()

  grid_cropped <- crop_background_spatial(background, region, crop = TRUE)
  bb <- sf::st_bbox(grid_cropped)
  region_bb <- sf::st_bbox(region)

  expect_true(bb["xmin"] >= region_bb["xmin"] - 1e-6)
  expect_true(bb["xmax"] <= region_bb["xmax"] + 1e-6)
  expect_true(bb["ymin"] >= region_bb["ymin"] - 1e-6)
  expect_true(bb["ymax"] <= region_bb["ymax"] + 1e-6)
})

test_that("crop_background_spatial() drops cells entirely outside region", {
  region <- make_test_region()
  background <- make_test_background()

  grid_cropped <- crop_background_spatial(background, region, crop = TRUE)
  full_grid <- background_to_grid_sf(background, region)

  expect_lt(nrow(grid_cropped), nrow(full_grid))
})

test_that("crop_background_spatial() works with a multi-polygon region via st_union", {
  region_multi <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(5, 0), c(5, 5), c(0, 5), c(0, 0)))),
    sf::st_polygon(list(rbind(c(6, 6), c(10, 6), c(10, 10), c(6, 10), c(6, 6)))),
    crs = 27700
  )
  background <- make_test_background()

  expect_no_error(
    crop_background_spatial(background, region_multi, crop = TRUE)
  )
})

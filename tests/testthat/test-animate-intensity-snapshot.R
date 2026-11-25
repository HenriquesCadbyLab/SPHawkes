make_snapshot_animation_fixture <- function() {
  x <- c(0, 1, 2)
  y <- c(0, 1, 2)

  background <- outer(
    x, y,
    function(x, y) x + 10 * y
  )

  list(
    days = c(10, 11),
    x = x,
    y = y,
    inside = rep(TRUE, length(x) * length(y)),
    snapshots = list(
      list(
        bg = background,
        tr = matrix(1, length(x), length(y))
      ),
      list(
        bg = 2 * background,
        tr = matrix(3, length(x), length(y))
      )
    )
  )
}

skip_without_animation_dependencies <- function() {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gifski")
}

test_that("animation writes a non-empty GIF and returns its filename invisibly", {
  skip_without_animation_dependencies()

  grid <- make_snapshot_animation_fixture()

  directory <- tempfile("spHawkes-animation-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  filename <- file.path(directory, "intensity.gif")

  result <- withVisible(
    animate_intensity_snapshots(
      grid,
      component = "total",
      filename = filename,
      fps = 2,
      width = 240,
      height = 200
    )
  )

  expect_identical(result$value, filename)
  expect_false(result$visible)
  expect_true(file.exists(filename))
  expect_gt(unname(file.info(filename)$size), 6)

  signature <- rawToChar(
    readBin(filename, what = "raw", n = 6L)
  )

  expect_true(signature %in% c("GIF87a", "GIF89a"))
})

test_that("animation accepts a subset of evaluated times", {
  skip_without_animation_dependencies()

  grid <- make_snapshot_animation_fixture()

  filename <- tempfile(fileext = ".gif")
  on.exit(unlink(filename), add = TRUE)

  expect_no_error(
    animate_intensity_snapshots(
      grid,
      component = "triggering",
      times = 11,
      filename = filename,
      width = 240,
      height = 200
    )
  )

  expect_true(file.exists(filename))
  expect_gt(unname(file.info(filename)$size), 6)
})

test_that("existing output is protected unless overwrite is enabled", {
  skip_without_animation_dependencies()

  grid <- make_snapshot_animation_fixture()

  filename <- tempfile(fileext = ".gif")
  on.exit(unlink(filename), add = TRUE)

  writeLines("existing output", filename)

  expect_error(
    animate_intensity_snapshots(grid, filename = filename),
    "output file exists"
  )

  expect_identical(readLines(filename), "existing output")

  expect_no_error(
    animate_intensity_snapshots(
      grid,
      times = 10,
      filename = filename,
      width = 240,
      height = 200,
      overwrite = TRUE
    )
  )

  signature <- rawToChar(
    readBin(filename, what = "raw", n = 6L)
  )

  expect_true(signature %in% c("GIF87a", "GIF89a"))
})

test_that("invalid animation arguments fail before writing output", {
  skip_without_animation_dependencies()

  grid <- make_snapshot_animation_fixture()

  filename <- tempfile(fileext = ".gif")
  on.exit(unlink(filename), add = TRUE)

  expect_error(
    animate_intensity_snapshots(list(), filename = filename),
    "must be an output"
  )
  expect_error(
    animate_intensity_snapshots(grid, times = 999, filename = filename),
    "must be available"
  )
  expect_error(
    animate_intensity_snapshots(
      grid, times = c(10, 10), filename = filename
    ),
    "distinct finite"
  )
  expect_error(
    animate_intensity_snapshots(grid, fps = 0, filename = filename),
    "'fps' must be positive",
    fixed = TRUE
  )
  expect_error(
    animate_intensity_snapshots(grid, width = 10.5, filename = filename),
    "positive integers"
  )
  expect_error(
    animate_intensity_snapshots(grid, overwrite = NA, filename = filename),
    "TRUE or FALSE"
  )
  expect_error(
    animate_intensity_snapshots(grid, filename = "intensity.png"),
    "GIF filename"
  )

  missing_directory <- tempfile("missing-animation-directory-")

  expect_error(
    animate_intensity_snapshots(
      grid,
      filename = file.path(missing_directory, "intensity.gif")
    ),
    "output directory does not exist"
  )

  expect_false(file.exists(filename))
})

test_that("invalid later snapshots are detected before export", {
  skip_without_animation_dependencies()

  grid <- make_snapshot_animation_fixture()
  grid$snapshots[[2]]$tr[1L] <- Inf

  filename <- tempfile(fileext = ".gif")
  on.exit(unlink(filename), add = TRUE)

  expect_error(
    animate_intensity_snapshots(
      grid,
      component = "triggering",
      filename = filename
    ),
    "finite and non-negative"
  )

  expect_false(file.exists(filename))
})

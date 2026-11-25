make_snapshot_plot_fixture <- function() {
  x <- c(0, 1, 2)
  y <- c(0, 1, 2, 3)

  background <- outer(
    x, y,
    function(x, y) x + 10 * y
  )

  inside <- rep(TRUE, length(x) * length(y))
  inside[1L] <- FALSE

  list(
    days = c(10, 11),
    x = x,
    y = y,
    inside = inside,
    snapshots = list(
      list(
        bg = background,
        tr = matrix(2, nrow = length(x), ncol = length(y))
      ),
      list(
        bg = 2 * background,
        tr = matrix(5, nrow = length(x), ncol = length(y))
      )
    )
  )
}

test_that("snapshot plots preserve coordinates and matrix orientation", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()

  p <- plot_intensity_snapshot(
    grid,
    time = 10,
    component = "background"
  )

  expect_s3_class(p, "ggplot")

  expected <- expand.grid(x = grid$x, y = grid$y)
  expected <- expected[grid$inside, , drop = FALSE]

  expect_equal(p$data$x, expected$x)
  expect_equal(p$data$y, expected$y)

  expect_equal(
    p$data$intensity,
    expected$x + 10 * expected$y
  )

  expect_equal(nrow(p$data), sum(grid$inside))
})

test_that("snapshot plots select model times rather than positions", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()

  p <- plot_intensity_snapshot(
    grid,
    time = 11,
    component = "background"
  )

  expect_equal(
    p$data$intensity,
    as.vector(grid$snapshots[[2]]$bg)[grid$inside]
  )

  expect_error(
    plot_intensity_snapshot(grid, time = 1),
    "not available"
  )
})

test_that("background, triggering and total components are consistent", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()

  background <- plot_intensity_snapshot(
    grid, time = 10, component = "background"
  )
  triggering <- plot_intensity_snapshot(
    grid, time = 10, component = "triggering"
  )
  total <- plot_intensity_snapshot(
    grid, time = 10, component = "total"
  )

  expect_equal(triggering$data$intensity, rep(2, sum(grid$inside)))

  expect_equal(
    total$data$intensity,
    background$data$intensity + triggering$data$intensity
  )

  expect_equal(
    total$scales$get_scales("fill")$limits,
    c(0, max(total$data$intensity))
  )
})

test_that("tile dimensions follow coordinate spacing or cell_size", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()
  grid$x <- c(0, 2, 4)
  grid$y <- c(0, 3, 6, 9)

  inferred <- ggplot2::ggplot_build(
    plot_intensity_snapshot(grid, time = 10)
  )$data[[1]]

  expect_equal(
    inferred$xmax - inferred$xmin,
    rep(2, nrow(inferred))
  )
  expect_equal(
    inferred$ymax - inferred$ymin,
    rep(3, nrow(inferred))
  )

  explicit <- ggplot2::ggplot_build(
    plot_intensity_snapshot(grid, time = 10, cell_size = 5)
  )$data[[1]]

  expect_equal(explicit$xmax - explicit$xmin, rep(5, nrow(explicit)))
  expect_equal(explicit$ymax - explicit$ymin, rep(5, nrow(explicit)))
})

test_that("zero intensities produce a valid colour scale", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()
  grid$snapshots[[1]]$bg[] <- 0
  grid$snapshots[[1]]$tr[] <- 0

  p <- plot_intensity_snapshot(grid, time = 10)

  expect_equal(p$data$intensity, rep(0, sum(grid$inside)))
  expect_equal(p$scales$get_scales("fill")$limits, c(0, 1))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("invalid snapshot inputs fail explicitly", {
  skip_if_not_installed("ggplot2")

  grid <- make_snapshot_plot_fixture()

  expect_error(
    plot_intensity_snapshot(list(), time = 10),
    "must be an output"
  )
  expect_error(
    plot_intensity_snapshot(grid, time = c(10, 11)),
    "single finite"
  )
  expect_error(
    plot_intensity_snapshot(grid, time = NA_real_),
    "single finite"
  )
  expect_error(
    plot_intensity_snapshot(grid, time = 10, component = "unknown")
  )
  expect_error(
    plot_intensity_snapshot(grid, time = 10, cell_size = 0),
    "positive finite"
  )

  bad <- grid
  bad$days <- c(10, 10)
  expect_error(
    plot_intensity_snapshot(bad, time = 10),
    "invalid times"
  )

  bad <- grid
  bad$inside[2L] <- NA
  expect_error(
    plot_intensity_snapshot(bad, time = 10),
    "invalid inside-region mask"
  )

  bad <- grid
  bad$inside[] <- FALSE
  expect_error(
    plot_intensity_snapshot(bad, time = 10),
    "No grid locations"
  )

  bad <- grid
  bad$snapshots[[1]]$bg <- matrix(1, 2, 2)
  expect_error(
    plot_intensity_snapshot(bad, time = 10),
    "invalid matrix dimensions"
  )

  bad <- grid
  bad$snapshots[[1]]$bg[2L] <- -1

  expect_error(
    plot_intensity_snapshot(
      bad,
      time = 10,
      component = "background"
    ),
    "finite and non-negative"
  )

  bad <- grid
  bad$x <- c(0, 1, 3)
  expect_error(
    plot_intensity_snapshot(bad, time = 10),
    "regularly spaced"
  )
})

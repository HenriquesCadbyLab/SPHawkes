make_integration_fixture <- function(times = c(0, 0.5, 1, 2)) {
  x <- c(0, 1)
  y <- c(0, 1, 2)

  spatial <- outer(
    x, y,
    function(x, y) 1 + x + 10 * y
  )

  inside <- rep(TRUE, length(x) * length(y))
  inside[1L] <- FALSE

  snapshots <- lapply(
    times,
    function(t) {
      list(
        bg = spatial * (1 + t),
        tr = spatial * (2 + 3 * t)
      )
    }
  )

  list(
    days = times,
    x = x,
    y = y,
    inside = inside,
    snapshots = snapshots
  )
}

expected_integration_surface <- function(grid) {
  spatial <- outer(
    grid$x, grid$y,
    function(x, y) 1 + x + 10 * y
  )

  spatial[
    !matrix(
      grid$inside,
      nrow = length(grid$x),
      ncol = length(grid$y)
    )
  ] <- 0

  spatial
}

test_that("integration returns a compatible grid with period metadata", {
  grid <- make_integration_fixture()

  result <- integrate_intensity_period(
    grid,
    start = 0,
    end = 2
  )

  expect_named(
    result,
    c(
      "days", "x", "y", "inside", "snapshots",
      "statistic", "period", "integration"
    )
  )

  expect_identical(result$x, grid$x)
  expect_identical(result$y, grid$y)
  expect_identical(result$inside, grid$inside)

  expect_equal(result$days, 0)
  expect_length(result$snapshots, 1L)
  expect_named(result$snapshots[[1]], c("bg", "tr"))

  expect_equal(
    dim(result$snapshots[[1]]$bg),
    c(length(grid$x), length(grid$y))
  )
  expect_equal(
    dim(result$snapshots[[1]]$tr),
    c(length(grid$x), length(grid$y))
  )

  expect_identical(result$statistic, "time_integral")
  expect_equal(result$period, c(start = 0, end = 2))
  expect_identical(result$integration$method, "trapezoid")
  expect_identical(result$integration$time_unit, "days")
  expect_equal(result$integration$times, grid$days)
})

test_that("linear temporal surfaces integrate exactly on irregular nodes", {
  grid <- make_integration_fixture(
    times = c(0, 0.25, 0.9, 1.4, 2)
  )

  result <- integrate_intensity_period(
    grid,
    start = 0,
    end = 2
  )

  spatial <- expected_integration_surface(grid)

  # Integral of 1 + t over [0, 2] is 4.
  expect_equal(result$snapshots[[1]]$bg, 4 * spatial)

  # Integral of 2 + 3t over [0, 2] is 10.
  expect_equal(result$snapshots[[1]]$tr, 10 * spatial)
})

test_that("constant surfaces are multiplied by period duration", {
  grid <- make_integration_fixture()

  background <- matrix(2, length(grid$x), length(grid$y))
  triggering <- matrix(3, length(grid$x), length(grid$y))

  grid$snapshots <- lapply(
    grid$days,
    function(t) list(bg = background, tr = triggering)
  )

  result <- integrate_intensity_period(
    grid,
    start = 0,
    end = 2
  )

  expected_bg <- matrix(4, length(grid$x), length(grid$y))
  expected_tr <- matrix(6, length(grid$x), length(grid$y))

  mask <- matrix(
    grid$inside,
    nrow = length(grid$x),
    ncol = length(grid$y)
  )

  expected_bg[!mask] <- 0
  expected_tr[!mask] <- 0

  expect_equal(result$snapshots[[1]]$bg, expected_bg)
  expect_equal(result$snapshots[[1]]$tr, expected_tr)
})

test_that("only the requested period contributes to integration", {
  grid <- make_integration_fixture()

  result <- integrate_intensity_period(
    grid,
    start = 0.5,
    end = 1
  )

  spatial <- expected_integration_surface(grid)

  expect_equal(result$snapshots[[1]]$bg, 0.875 * spatial)
  expect_equal(result$snapshots[[1]]$tr, 2.125 * spatial)
  expect_equal(result$integration$times, c(0.5, 1))
  expect_equal(result$days, 0.5)
})

test_that("unsorted snapshots are reordered together with their times", {
  grid <- make_integration_fixture()

  reference <- integrate_intensity_period(grid, start = 0, end = 2)

  permutation <- c(4, 2, 1, 3)
  grid$days <- grid$days[permutation]
  grid$snapshots <- grid$snapshots[permutation]

  result <- integrate_intensity_period(grid, start = 0, end = 2)

  expect_equal(result, reference)
})

test_that("integration is additive across adjacent periods", {
  grid <- make_integration_fixture()

  whole <- integrate_intensity_period(grid, start = 0, end = 2)
  first <- integrate_intensity_period(grid, start = 0, end = 1)
  second <- integrate_intensity_period(grid, start = 1, end = 2)

  for (component in c("bg", "tr")) {
    expect_equal(
      whole$snapshots[[1]][[component]],
      first$snapshots[[1]][[component]] +
        second$snapshots[[1]][[component]]
    )
  }
})

test_that("outside-region locations are zeroed without modifying input", {
  grid <- make_integration_fixture()

  for (i in seq_along(grid$snapshots)) {
    grid$snapshots[[i]]$bg[1L] <- NA_real_
    grid$snapshots[[i]]$tr[1L] <- Inf
  }

  original <- grid
  result <- integrate_intensity_period(grid, start = 0, end = 2)

  expect_equal(result$snapshots[[1]]$bg[1L], 0)
  expect_equal(result$snapshots[[1]]$tr[1L], 0)
  expect_identical(grid, original)

  expect_true(all(is.finite(result$snapshots[[1]]$bg)))
  expect_true(all(is.finite(result$snapshots[[1]]$tr)))
})

test_that("two boundary snapshots are sufficient", {
  grid <- make_integration_fixture(times = c(0, 2))

  result <- integrate_intensity_period(grid, start = 0, end = 2)
  spatial <- expected_integration_surface(grid)

  expect_equal(result$snapshots[[1]]$bg, 4 * spatial)
  expect_equal(result$snapshots[[1]]$tr, 10 * spatial)
})

test_that("integrated surfaces can be passed to the snapshot plotter", {
  skip_if_not_installed("ggplot2")

  grid <- make_integration_fixture()
  result <- integrate_intensity_period(grid, start = 0, end = 2)

  p <- plot_intensity_snapshot(
    grid = result,
    time = 0,
    component = "total"
  )

  expected_total <- result$snapshots[[1]]$bg +
    result$snapshots[[1]]$tr

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data$intensity,
    as.vector(expected_total)[result$inside]
  )
})

test_that("invalid periods and unavailable boundaries are rejected", {
  grid <- make_integration_fixture()

  expect_error(
    integrate_intensity_period(grid, start = 1, end = 1),
    "'end' > 'start'",
    fixed = TRUE
  )
  expect_error(
    integrate_intensity_period(grid, start = 2, end = 0),
    "'end' > 'start'",
    fixed = TRUE
  )
  expect_error(
    integrate_intensity_period(grid, start = NA_real_, end = 2),
    "finite times"
  )
  expect_error(
    integrate_intensity_period(grid, start = c(0, 1), end = 2),
    "finite times"
  )
  expect_error(
    integrate_intensity_period(grid, start = 0, end = Inf),
    "finite times"
  )
  expect_error(
    integrate_intensity_period(grid, start = 0, end = 1.5),
    "Both period boundaries"
  )
})

test_that("invalid grids and repeated integration are rejected", {
  grid <- make_integration_fixture()

  expect_error(
    integrate_intensity_period(list(), start = 0, end = 2),
    "must be an output"
  )

  bad <- grid
  bad$days[2L] <- bad$days[1L]
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "distinct valid times"
  )

  bad <- grid
  bad$snapshots <- bad$snapshots[-1L]
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "distinct valid times"
  )

  bad <- grid
  bad$x <- rev(bad$x)
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "strictly increasing"
  )

  bad <- grid
  bad$inside[2L] <- NA
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "inside-region mask"
  )

  bad <- grid
  bad$inside[] <- FALSE
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "inside-region mask"
  )

  integrated <- integrate_intensity_period(grid, start = 0, end = 2)

  expect_error(
    integrate_intensity_period(integrated, start = 0, end = 2),
    "Only instantaneous snapshots"
  )
})

test_that("invalid component matrices and inside-region values are rejected", {
  grid <- make_integration_fixture()

  bad <- grid
  bad$snapshots[[2L]]$bg <- matrix(1, 2, 2)
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "invalid intensity matrices"
  )

  bad <- grid
  bad$snapshots[[2L]]$tr <- NULL
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "invalid intensity matrices"
  )

  # Location 2 is inside the region.
  bad <- grid
  bad$snapshots[[2L]]$bg[2L] <- -1
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "finite and non-negative"
  )

  bad <- grid
  bad$snapshots[[2L]]$tr[2L] <- Inf
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "finite and non-negative"
  )

  bad <- grid
  bad$snapshots[[2L]]$bg[2L] <- NA_real_
  expect_error(
    integrate_intensity_period(bad, start = 0, end = 2),
    "finite and non-negative"
  )
})

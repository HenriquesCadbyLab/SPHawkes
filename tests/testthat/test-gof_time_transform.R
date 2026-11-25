fit <- make_gof_fit()

test_that("gof_time_transform() validates its inputs", {
  expect_error(
    gof_time_transform(fit, dt = -1),
    "`dt` must be a single positive finite number"
  )

  expect_error(
    gof_time_transform(fit, dt = 0),
    "`dt` must be a single positive finite number"
  )

  expect_error(
    gof_time_transform(fit, level = 1),
    "`level` must be a single number strictly between 0 and 1"
  )

  expect_error(
    gof_time_transform(fit, level = 0),
    "`level` must be a single number strictly between 0 and 1"
  )
})

test_that("gof_time_transform() returns the expected structure", {
  gof <- gof_time_transform(fit, dt = 5, level = 0.95)

  expect_s3_class(gof, "semihawkes_gof_time")
  expect_type(gof$grid, "list")
  expect_type(gof$events, "list")
  expect_type(gof$bands, "list")

  expect_named(
    gof$grid,
    c(
      "time",
      "transformed_time",
      "observed_cumulative",
      "expected_cumulative"
    )
  )

  expect_named(
    gof$events,
    c("days", "transformed_time")
  )

  expect_named(
    gof$bands,
    c("count", "lower", "upper")
  )

  expect_identical(gof$level, 0.95)
  expect_identical(gof$dt, 5)
})

test_that("cumulative counts and transformed times are coherent", {
  gof <- gof_time_transform(fit, dt = 5)

  n_events <- nrow(fit$events$events)

  expect_identical(
    max(gof$grid$observed_cumulative),
    as.numeric(n_events)
  )

  expect_identical(
    max(gof$bands$count),
    n_events
  )

  expect_true(all(diff(gof$grid$transformed_time) >= 0))
  expect_true(all(diff(gof$events$transformed_time) >= 0))

  expect_equal(
    max(gof$grid$expected_cumulative),
    fit$lambda_at_all,
    tolerance = 0.05
  )
})

test_that("KS bands are valid and ordered", {
  gof <- gof_time_transform(fit, dt = 5, level = 0.95)

  n_events <- nrow(fit$events$events)

  expect_identical(gof$bands$count, 0:n_events)

  expect_true(all(gof$bands$lower >= 0))
  expect_true(all(gof$bands$upper <= n_events))
  expect_true(all(gof$bands$lower <= gof$bands$upper))

  expect_true(all(diff(gof$bands$lower) >= 0))
  expect_true(all(diff(gof$bands$upper) >= 0))
})

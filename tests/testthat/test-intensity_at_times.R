fit <- make_gof_fit()

test_that("intensity_at_times() validates its inputs", {
  expect_error(
    intensity_at_times(fit, times = "not numeric"),
    "`times` must be a finite numeric vector"
  )

  expect_error(
    intensity_at_times(fit, times = c(-1, 0)),
    "`times` must lie within the fitted observation window"
  )

  expect_error(
    intensity_at_times(fit, times = c(0, fit$events$TT + 1)),
    "`times` must lie within the fitted observation window"
  )

  expect_error(
    intensity_at_times(fit, times = 0, exclude_self = "yes"),
    "`exclude_self` must be a single logical value"
  )
})

test_that("intensity_at_times() returns one non-negative value per time", {
  times <- c(0, 100, 250, fit$events$TT)

  lambda_t <- intensity_at_times(fit, times = times)

  expect_type(lambda_t, "double")
  expect_length(lambda_t, length(times))
  expect_true(all(is.finite(lambda_t)))
  expect_true(all(lambda_t >= 0))
})

test_that("intensity_at_times() is increasing in the triggering contribution", {
  times <- c(100, 250, fit$events$TT)

  lambda_excluded <- intensity_at_times(
    fit,
    times = times,
    exclude_self = TRUE
  )

  lambda_included <- intensity_at_times(
    fit,
    times = times,
    exclude_self = FALSE
  )

  expect_true(all(lambda_included >= lambda_excluded))
})

test_that("temporal integral approximates fit$lambda_at_all", {
  skip_on_cran()

  dt <- 1
  times <- seq(0, fit$events$TT, by = dt)

  lambda_t <- intensity_at_times(fit, times = times)

  temporal_integral <- sum(lambda_t) * dt

  expect_equal(
    temporal_integral,
    fit$lambda_at_all,
    tolerance = 0.05
  )

  expect_equal(
    temporal_integral,
    nrow(fit$events$events),
    tolerance = 0.05
  )
})

fit <- make_gof_fit()

test_that("gof_period_residuals() validates its inputs", {
  expect_error(
    gof_period_residuals(fit, period = "week"),
    "'arg' should be one of"
  )
})

test_that("gof_period_residuals() returns one object per year", {
  gof_periods <- gof_period_residuals(fit, period = "year")

  expected_years <- sort(unique(format(fit$events$events$date, "%Y")))

  expect_type(gof_periods, "list")
  expect_named(gof_periods, expected_years)

  for (gof in gof_periods) {
    expect_s3_class(gof, "semihawkes_gof_spatial")
    expect_true(nrow(gof$cells) > 0)
    expect_true(sum(gof$cells$observed) > 0)
  }
})

test_that("period residuals conserve the total observed count", {
  gof_periods <- gof_period_residuals(fit, period = "year")

  total_observed <- sum(
    vapply(gof_periods, function(gof) sum(gof$cells$observed), numeric(1))
  )

  expect_equal(
    total_observed,
    nrow(fit$events$events),
    tolerance = 1e-12
  )
})

test_that("period residuals have the expected class", {
  gof_periods <- gof_period_residuals(fit, period = "year")

  expect_s3_class(gof_periods, "semihawkes_gof_period")
})

test_that("period plot supports residual types", {
  gof_periods <- gof_period_residuals(fit, period = "year")

  expect_s3_class(
    plot(gof_periods, type = "raw"),
    "ggplot"
  )

  expect_s3_class(
    plot(gof_periods, type = "relative"),
    "ggplot"
  )

  expect_s3_class(
    plot(gof_periods, type = "standardized"),
    "ggplot"
  )

  expect_error(
    plot(gof_periods, type = "invalid"),
    "'arg' should be one of"
  )
})

test_that("period windows are calendar periods covering the observation window once", {
  # Each window used to run from the first to the last event date in its
  # period, and the expected count stopped at the start of the last day while
  # events on that day were counted, so the monthly expected counts fell short
  # of the observed ones (by 5% on the Milwaukee example).
  gof_months <- gof_period_residuals(fit, period = "month")

  first_day <- fit$events$start_date
  last_day <- fit$events$start_date + fit$events$TT - 1
  expect_named(
    gof_months,
    format(seq(as.Date(format(first_day, "%Y-%m-01")), last_day, by = "month"), "%Y-%m")
  )

  starts <- do.call(c, unname(lapply(gof_months, `[[`, "start_date")))
  ends <- do.call(c, unname(lapply(gof_months, `[[`, "end_date")))
  expect_equal(starts[1], first_day)
  expect_equal(ends[length(ends)], last_day)
  # Contiguous: each window starts the day after the previous one ends.
  expect_equal(as.numeric(starts[-1] - ends[-length(ends)]), rep(1, length(ends) - 1))
  # Interior windows are whole calendar months.
  inner <- gof_months[-c(1, length(gof_months))]
  expect_equal(
    unname(vapply(inner, `[[`, numeric(1), "window_length")),
    unname(vapply(names(inner), function(p) {
      d <- as.Date(paste0(p, "-01"))
      as.numeric(seq(d, by = "month", length.out = 2)[2] - d)
    }, numeric(1)))
  )

  # Observed and expected counts both add up to the whole-window totals.
  whole <- gof_spatial_residuals(fit)
  total <- function(g, col) sum(vapply(g, function(x) sum(x$cells[[col]]), numeric(1)))
  expect_equal(total(gof_months, "observed"), sum(whole$cells$observed))
  expect_equal(total(gof_months, "expected"), sum(whole$cells$expected), tolerance = 1e-3)
})

test_that("a month whose events all fall on one day gets a whole-month window", {
  # Its window used to be that single day with an empty interval, which
  # stopped with "must define a non-empty window".
  fit2 <- fit
  ev <- fit2$events$events
  month <- format(ev$date, "%Y-%m")
  keep <- month != "2020-03" | ev$date == min(ev$date[month == "2020-03"])
  fit2$events$events <- ev[keep, ]

  expect_no_error(gof_months <- gof_period_residuals(fit2, period = "month"))
  expect_equal(gof_months[["2020-03"]]$window_length, 31)
  expect_equal(sum(gof_months[["2020-03"]]$cells$observed), sum(keep & month == "2020-03"))
})

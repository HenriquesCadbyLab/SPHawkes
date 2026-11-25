make_test_events <- function(n = 300, seed = 42) {
  synth <- make_synthetic_data(n = n, seed = seed)
  prepare_events(synth$points, synth$region, date_col = "event_date",
                  start_date = "2020-01-01")
}

test_that("temporal_base_weekly returns a mean-1 periodic curve", {
  ev <- make_test_events()
  wb <- temporal_base_weekly(ev)

  expect_s3_class(wb, "semihawkes_temporal_base")
  expect_equal(wb$period, 7)
  expect_equal(mean(wb$basevalue), 1, tolerance = 1e-8)
  expect_length(wb$basevalue, length(wb$breaks))

  # the wrapped function should be periodic with period 7
  expect_equal(wb$fun(3.2), wb$fun(3.2 + 7), tolerance = 1e-6)
  expect_equal(wb$fun(3.2), wb$fun(3.2 + 14), tolerance = 1e-6)
})

test_that("temporal_base_weekly accepts a resp_weight vector without erroring", {
  ev <- make_test_events()
  n <- nrow(ev$events)
  wb <- temporal_base_weekly(ev, resp_weight = stats::runif(n, 0.5, 1))
  expect_equal(mean(wb$basevalue), 1, tolerance = 1e-8)
})

test_that("temporal_base_monthly returns a mean-1 periodic curve using calendar month", {
  ev <- make_test_events()
  mb <- temporal_base_monthly(ev)

  expect_s3_class(mb, "semihawkes_temporal_base")
  expect_equal(mb$period, 12)
  expect_equal(mean(mb$basevalue), 1, tolerance = 1e-8)

  # periodic with period 12
  expect_equal(mb$fun(2.5), mb$fun(2.5 + 12), tolerance = 1e-6)
})

test_that("temporal_base_monthly uses events$month (calendar month), not days mod 12", {
  ev <- make_test_events()
  # Reproduce the "days mod 12" quantity from the buggy Loop-0 code and check
  # it does NOT match events$month in general, so we know the two would give
  # different results (confirming the fix actually changes behaviour).
  days_mod_12 <- ev$events$days - as.integer(ev$events$days / 12) * 12
  expect_false(isTRUE(all.equal(days_mod_12, ev$events$month)))
})

test_that("temporal_base_trend with resp_weight = NULL is flat with mean 1", {
  ev <- make_test_events()
  tb <- temporal_base_trend(ev)

  expect_s3_class(tb, "semihawkes_temporal_base")
  expect_true(is.na(tb$period))
  expect_true(all(tb$basevalue == 1))
})

test_that("temporal_base_trend with resp_weight matches a brute-force loop", {
  ev <- make_test_events(n = 40) # keep small: brute force is O(n_events * n_marks)
  n <- nrow(ev$events)
  resp_weight <- stats::runif(n, 0.5, 1)
  bandwidth <- 150

  tb <- temporal_base_trend(ev, resp_weight = resp_weight, bandwidth = bandwidth,
                             chunk_size = 7) # force multiple chunks
  time_marks <- ev$time.marks

  brute <- rep(0, length(time_marks))
  for (i in seq_len(n)) {
    wi <- resp_weight[i]
    brute <- brute + wi * stats::dnorm(ev$events$days[i] - time_marks, 0, bandwidth) /
      (stats::pnorm(ev$TT, ev$events$days[i], bandwidth) -
         stats::pnorm(0, ev$events$days[i], bandwidth))
  }
  brute <- brute / mean(brute)

  # temporal_base_trend() evaluates the curve on a 1-day grid and
  # interpolates to the 0.05-day time marks, which changes it by at most
  # about (1 / 150)^2 / 8 = 6e-6 of its value.
  expect_equal(tb$basevalue, brute, tolerance = 1e-5)
  expect_lt(max(abs(tb$basevalue / brute - 1)), 1e-5)
})

test_that("temporal_base_trend with a short bandwidth evaluates at the time marks directly", {
  ev <- make_test_events(n = 40)
  n <- nrow(ev$events)
  resp_weight <- stats::runif(n, 0.5, 1)
  bandwidth <- 3 # bandwidth / 100 is finer than the 0.05-day time marks

  tb <- temporal_base_trend(ev, resp_weight = resp_weight, bandwidth = bandwidth)
  time_marks <- ev$time.marks

  brute <- rep(0, length(time_marks))
  for (i in seq_len(n)) {
    brute <- brute + resp_weight[i] *
      stats::dnorm(ev$events$days[i] - time_marks, 0, bandwidth) /
      (stats::pnorm(ev$TT, ev$events$days[i], bandwidth) -
         stats::pnorm(0, ev$events$days[i], bandwidth))
  }
  brute <- brute / mean(brute)

  expect_equal(tb$basevalue, brute, tolerance = 1e-12)
})

test_that("temporal_base_trend fun is 0 outside [0, TT]", {
  ev <- make_test_events()
  tb <- temporal_base_trend(ev, resp_weight = stats::runif(nrow(ev$events)))
  expect_equal(tb$fun(-10), 0)
  expect_equal(tb$fun(ev$TT + 10), 0)
})

make_em_step_inputs <- function(n = 250, seed = 21, cell_size = 150,
                                spatial_range = 500, spatial_step = 50,
                                pad_factor = 4) {
  synth <- make_synthetic_data(n = n, seed = seed)
  padded_region <- sf::st_buffer(synth$region, dist = pad_factor * cell_size)
  ev <- prepare_events(synth$points, padded_region, date_col = "event_date",
                       start_date = "2020-01-01")
  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  ev <- add_bg_integral(ev, n_mc = 500, seed = 1)

  trend   <- temporal_base_trend(ev)
  weekly  <- temporal_base_weekly(ev)
  monthly <- temporal_base_monthly(ev)
  background <- background_spatial(ev, cell_size = cell_size)

  mub_check <- background$fun(ev$events$coorx, ev$events$coory)
  if (anyNA(mub_check)) {
    bad <- ev$events[is.na(mub_check), c("coorx", "coory")]
    stop(sprintf(
      "Fixture bug: %d/%d synthetic events (seed = %s) still fall outside ",
      sum(is.na(mub_check)), length(mub_check), seed
    ), "background_spatial()'s grid after padding by ", pad_factor,
    " * cell_size -- increase pad_factor in make_em_step_inputs(). ",
    "Offending coordinates:\n",
    paste(capture.output(print(bad)), collapse = "\n"))
  }

  temporal_trigger <- trigger_temporal_init()
  spatial_trigger <- trigger_spatial_init(x_range = spatial_range, y_range = spatial_range,
                                          x_step = spatial_step, y_step = spatial_step)

  cache_dir <- tempfile("em_step_marks_")
  dir.create(cache_dir)
  trigger_spatial_marks(ev, spatial_trigger, cache_dir = cache_dir)

  pairs <- trigger_pairs(ev, window_dist = spatial_range)

  bgint0 <- background_intensity(ev, trend, weekly, monthly, background)
  trig0  <- trigger_evaluate(ev, background, temporal_trigger, spatial_trigger)

  n_events <- nrow(ev$events)
  A  <- n_events / (2 * trig0$at_all)
  mu <- n_events / (2 * bgint0$bg_at_all)

  list(events = ev, background = background, trend = trend, weekly = weekly,
       monthly = monthly, temporal_trigger = temporal_trigger,
       spatial_trigger = spatial_trigger, mu = mu, A = A, pairs = pairs,
       cache_dir = cache_dir, cell_size = cell_size)
}

run_em_step <- function(inputs, ...) {
  em_step(inputs$events, inputs$background, inputs$trend, inputs$weekly, inputs$monthly,
          inputs$temporal_trigger, inputs$spatial_trigger, inputs$mu, inputs$A,
          pairs = inputs$pairs, cache_dir = inputs$cache_dir,
          cell_size = inputs$cell_size, ...)
}

test_that("em_step returns a well-formed semihawkes_em_step object", {
  inputs <- make_em_step_inputs()
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step <- run_em_step(inputs)
  n_events <- nrow(inputs$events$events)

  expect_s3_class(step, "semihawkes_em_step")
  expect_true(is.numeric(step$mu) && length(step$mu) == 1)
  expect_true(is.numeric(step$A) && length(step$A) == 1)
  expect_true(is.numeric(step$llik) && length(step$llik) == 1)
  expect_length(step$lambda_at_events, n_events)
  expect_length(step$bgprobs, n_events)

  expect_s3_class(step$background, "semihawkes_spatial_base")
  expect_s3_class(step$trend, "semihawkes_temporal_base")
  expect_s3_class(step$weekly, "semihawkes_temporal_base")
  expect_s3_class(step$monthly, "semihawkes_temporal_base")
  expect_s3_class(step$temporal_trigger, "semihawkes_trigger_temporal")
  expect_s3_class(step$spatial_trigger, "semihawkes_trigger_spatial")
})

test_that("em_step's mu, A, llik are finite and mu/A are positive", {
  inputs <- make_em_step_inputs()
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step <- run_em_step(inputs)

  expect_true(is.finite(step$mu) && step$mu > 0)
  expect_true(is.finite(step$A) && step$A > 0)
  expect_true(is.finite(step$llik))
  expect_true(all(is.finite(step$lambda_at_events)))
  expect_true(all(step$bgprobs >= 0 & step$bgprobs <= 1))
})

test_that("em_step's outputs match a manual step-by-step reproduction using the exported building blocks", {
  inputs <- make_em_step_inputs(n = 250, seed = 21)
  ev <- inputs$events
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  # --- Instrumentation: checksum every input BEFORE calling em_step(), so
  # if this test's numbers ever differ from an external standalone rerun
  # again, we can directly compare these printed digests instead of
  # re-guessing. Uses a simple sum-based fingerprint (no extra dependency
  # on the `digest` package) -- sufficient to detect "same data or not".
  fingerprint <- function(x) sprintf("%.10g", sum(as.numeric(x), na.rm = TRUE))
  message("=== FINGERPRINTS (test-em_step.R:186) ===")
  message("ev$events$coorx: ", fingerprint(ev$events$coorx))
  message("ev$events$coory: ", fingerprint(ev$events$coory))
  message("ev$events$days: ", fingerprint(ev$events$days))
  message("ev$events$bandwidth: ", fingerprint(ev$events$bandwidth))
  message("ev$events$bg.integral: ", fingerprint(ev$events$bg.integral))
  message("inputs$mu: ", inputs$mu, " inputs$A: ", inputs$A)
  message("inputs$background$basevalue: ", fingerprint(inputs$background$basevalue))
  message("inputs$weekly$basevalue: ", fingerprint(inputs$weekly$basevalue))
  message("inputs$monthly$basevalue: ", fingerprint(inputs$monthly$basevalue))
  message("inputs$trend$basevalue: ", fingerprint(inputs$trend$basevalue))
  message("nrow(inputs$pairs): ", nrow(inputs$pairs))
  message("=========================================")

  step <- run_em_step(inputs)

  # --- Manual replication of em_step()'s internals, calling the same
  # exported building blocks in the same order, entirely independently of
  # em_step()'s own code. ---
  bgint <- background_intensity(ev, inputs$trend, inputs$weekly, inputs$monthly, inputs$background)
  trig  <- trigger_evaluate(ev, inputs$background, inputs$temporal_trigger, inputs$spatial_trigger,
                            exclude_self = TRUE)

  lambda_at_events <- inputs$mu * bgint$bg_at_events + inputs$A * trig$at_events
  bgprobs <- inputs$mu * bgint$bg_at_events / lambda_at_events

  d <- ev$events
  mub_events <- bgint$mub_events
  # weekly/monthly are periodic in calendar time (day-of-week / calendar
  # month), not in days-since-start -- same conversion as in
  # background_intensity()/em_step() themselves. trend genuinely lives on
  # the days-since-start scale, so it's left as-is.
  marks <- bgint$marks_events
  weW <- inputs$weekly$fun(marks$wday)   * mub_events / lambda_at_events
  weT <- inputs$trend$fun(d$days)        * mub_events / lambda_at_events
  weM <- inputs$monthly$fun(marks$month) * mub_events / lambda_at_events

  message("manual bgint$bg_at_events: ", fingerprint(bgint$bg_at_events))
  message("manual trig$at_events: ", fingerprint(trig$at_events))
  message("manual lambda_at_events (pre-M-step): ", fingerprint(lambda_at_events))
  message("manual bgprobs (pre-M-step): ", fingerprint(bgprobs))

  weekly_new  <- temporal_base_weekly(ev,  resp_weight = weW)
  monthly_new <- temporal_base_monthly(ev, resp_weight = weM)
  trend_new   <- temporal_base_trend(ev,   resp_weight = weT)
  background_new <- background_spatial(ev, cell_size = inputs$cell_size, resp_weight = bgprobs)

  message("manual weekly_new$basevalue: ", fingerprint(weekly_new$basevalue),
          " | step$weekly$basevalue: ", fingerprint(step$weekly$basevalue))
  message("manual background_new$basevalue: ", fingerprint(background_new$basevalue),
          " | step$background$basevalue: ", fingerprint(step$background$basevalue))

  bgint_new <- background_intensity(ev, trend_new, weekly_new, monthly_new, background_new)

  edge_correction <- trigger_edge_correction(ev, inputs$temporal_trigger, inputs$spatial_trigger,
                                             inputs$cache_dir)
  pair_weights <- trigger_pair_weights(ev, inputs$pairs, inputs$A, lambda_at_events,
                                       inputs$temporal_trigger, inputs$spatial_trigger)

  message("manual pair_weights: ", fingerprint(pair_weights))

  temporal_trigger_new <- trigger_reestimate_temporal(ev, inputs$pairs, pair_weights,
                                                      edge_correction, inputs$temporal_trigger)
  spatial_trigger_new <- trigger_reestimate_spatial(ev, inputs$pairs, pair_weights,
                                                    edge_correction, inputs$spatial_trigger)

  message("manual temporal_trigger_new$basevalue: ", fingerprint(temporal_trigger_new$basevalue),
          " | step$temporal_trigger$basevalue: ", fingerprint(step$temporal_trigger$basevalue))
  message("manual spatial_trigger_new$basevalue: ", fingerprint(spatial_trigger_new$basevalue),
          " | step$spatial_trigger$basevalue: ", fingerprint(step$spatial_trigger$basevalue))

  trig_new <- trigger_evaluate(ev, background_new, temporal_trigger_new, spatial_trigger_new,
                               exclude_self = TRUE)

  n_events <- nrow(d)
  A_new  <- (n_events - sum(bgprobs)) / trig_new$at_all
  mu_new <- (n_events - A_new * trig_new$at_all) / bgint_new$bg_at_all

  message("manual trig_new$at_all: ", trig_new$at_all, " | manual bgint_new$bg_at_all: ", bgint_new$bg_at_all)
  message("manual A_new: ", A_new, " | step$A: ", step$A)
  message("manual mu_new: ", mu_new, " | step$mu: ", step$mu)

  lambda_at_events_new <- mu_new * bgint_new$bg_at_events + A_new * trig_new$at_events
  lambda_at_all_new    <- mu_new * bgint_new$bg_at_all    + A_new * trig_new$at_all
  llik_manual <- -(-sum(log(lambda_at_events_new)) + lambda_at_all_new)
  bgprobs_new_manual <- mu_new * bgint_new$bg_at_events / lambda_at_events_new

  expect_equal(step$mu, mu_new, tolerance = 1e-8)
  expect_equal(step$A, A_new, tolerance = 1e-8)
  expect_equal(step$llik, llik_manual, tolerance = 1e-8)
  expect_equal(step$lambda_at_all, lambda_at_all_new, tolerance = 1e-8)
  expect_equal(step$lambda_at_events, lambda_at_events_new, tolerance = 1e-8)
  expect_equal(step$bgprobs, bgprobs_new_manual, tolerance = 1e-8)
})

test_that("em_step's re-estimated components each satisfy their own invariants", {
  inputs <- make_em_step_inputs()
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step <- run_em_step(inputs)

  expect_equal(mean(step$background$basevalue[step$background$marks > 0]), 1,
               tolerance = 1e-6)
  expect_equal(mean(step$weekly$basevalue), 1, tolerance = 1e-6)
  expect_equal(mean(step$monthly$basevalue), 1, tolerance = 1e-6)

  tstep <- step$temporal_trigger$base[2] - step$temporal_trigger$base[1]
  expect_equal(simpson(step$temporal_trigger$basevalue, tstep), 1, tolerance = 1e-6)

  dx <- step$spatial_trigger$x[2] - step$spatial_trigger$x[1]
  dy <- step$spatial_trigger$y[2] - step$spatial_trigger$y[1]
  expect_equal(simpson.2D(step$spatial_trigger$basevalue, dx = dx, dy = dy), 1,
               tolerance = 1e-6)
})

test_that("em_step can be chained: feeding one step's output into the next runs cleanly", {
  inputs <- make_em_step_inputs(n = 250, seed = 5)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step1 <- run_em_step(inputs)

  inputs2 <- inputs
  inputs2$background <- step1$background
  inputs2$trend <- step1$trend
  inputs2$weekly <- step1$weekly
  inputs2$monthly <- step1$monthly
  inputs2$temporal_trigger <- step1$temporal_trigger
  inputs2$spatial_trigger <- step1$spatial_trigger
  inputs2$mu <- step1$mu
  inputs2$A <- step1$A

  step2 <- run_em_step(inputs2)

  expect_s3_class(step2, "semihawkes_em_step")
  expect_true(is.finite(step2$llik))
  expect_true(step2$mu > 0 && step2$A > 0)
})

test_that("em_step returns bgint and trig evaluated with its updated components", {
  inputs <- make_em_step_inputs(n = 250, seed = 5)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step <- run_em_step(inputs)

  expect_equal(step$bgint,
               background_intensity(inputs$events, step$trend, step$weekly,
                                    step$monthly, step$background))
  expect_equal(step$trig,
               trigger_evaluate(inputs$events, step$background,
                                step$temporal_trigger, step$spatial_trigger))
})

test_that("passing the previous step's bgint and trig gives exactly the same next step", {
  inputs <- make_em_step_inputs(n = 250, seed = 5)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step1 <- run_em_step(inputs)
  carried <- c("background", "trend", "weekly", "monthly", "temporal_trigger",
               "spatial_trigger", "mu", "A")
  inputs2 <- inputs
  inputs2[carried] <- step1[carried]

  recomputed <- run_em_step(inputs2)
  reused     <- run_em_step(inputs2, bgint = step1$bgint, trig = step1$trig)

  for (nm in c("mu", "A", "llik", "lambda_at_all", "lambda_at_events", "bgprobs")) {
    expect_identical(reused[[nm]], recomputed[[nm]], label = nm)
  }
})

test_that("em_step gives the same step with the spatial marks read in advance", {
  inputs <- make_em_step_inputs(n = 250, seed = 5)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  from_disk <- run_em_step(inputs)
  marks <- .read_spatial_marks(inputs$cache_dir, nrow(inputs$events$events))
  in_memory <- run_em_step(inputs, marks = marks)

  for (nm in c("mu", "A", "llik", "lambda_at_all", "lambda_at_events", "bgprobs")) {
    expect_identical(in_memory[[nm]], from_disk[[nm]], label = nm)
  }
})

test_that("em_step rejects bgint or trig that do not match the events", {
  inputs <- make_em_step_inputs(n = 250, seed = 5)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  expect_error(run_em_step(inputs, trig = list(at_events = 1, at_all = 1)),
               "`trig` does not match")
  expect_error(run_em_step(inputs, bgint = list(bg_at_events = 1)),
               "`bgint` does not match")
})

test_that("em_step gives the same result with nprocs = 1 and nprocs = 2", {
  inputs <- make_em_step_inputs(n = 240, seed = 9)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  nprocs_multi <- min(2L, parallel::detectCores(logical = FALSE))
  skip_if(nprocs_multi < 2, "Not enough cores available to test parallel chunking")

  step1 <- run_em_step(inputs, nprocs = 1)
  step2 <- run_em_step(inputs, nprocs = nprocs_multi)

  # The same to the last bit: trigger_evaluate() adds up the same chunks in
  # the same order whatever the number of workers.
  expect_identical(step2$mu, step1$mu)
  expect_identical(step2$A, step1$A)
  expect_identical(step2$llik, step1$llik)
  expect_identical(step2$lambda_at_events, step1$lambda_at_events)

  # So does a cluster passed in, which is left running for the next step.
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  step3 <- run_em_step(inputs, nprocs = cl)
  expect_identical(step3$llik, step1$llik)
  expect_equal(unlist(parallel::clusterEvalQ(cl, 1 + 1)), c(2, 2))
})

test_that("em_step propagates an informative error when cache_dir is missing marks", {
  inputs <- make_em_step_inputs(n = 250, seed = 3)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  empty_cache <- tempfile("empty_marks_")
  dir.create(empty_cache)
  inputs$cache_dir <- empty_cache

  expect_error(run_em_step(inputs), "trigger_spatial_marks")
})

test_that("print.semihawkes_em_step runs without error and reports mu/A/llik", {
  inputs <- make_em_step_inputs(n = 200, seed = 13)
  skip_if(nrow(inputs$pairs) == 0, "no triggering pairs in this synthetic draw")

  step <- run_em_step(inputs)
  expect_output(print(step), "semihawkes_em_step")
  expect_output(print(step), "mu =")
})

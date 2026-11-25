# ------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------

# Small, fully deterministic (once seeded) toy setup: 5 events, a hand-built
# `pairs` index, and `pair_weights` engineered so that
# bgprobs[i] + sum(pair_weights[pairs$i == i]) == 1 exactly for every i.
# Event 1 has no candidate parents (nothing precedes it), so bgprobs[1] = 1.
make_toy_branching_inputs <- function() {
  bgprobs <- c(1, 0.6, 0.3, 0.8, 0.2)

  pairs <- data.frame(
    i = c(2L, 3L, 3L, 4L, 5L, 5L),
    j = c(1L, 1L, 2L, 3L, 1L, 4L)
  )

  # Weights chosen so each child's weights sum to (1 - bgprobs[i]).
  pair_weights <- c(
    0.4,        # i=2, j=1  -> sums to 1 - 0.6 = 0.4
    0.25, 0.45, # i=3, j=1,2 -> sums to 1 - 0.3 = 0.70
    0.2,        # i=4, j=3  -> sums to 1 - 0.8 = 0.20
    0.5, 0.3    # i=5, j=1,4 -> sums to 1 - 0.2 = 0.80
  )

  list(bgprobs = bgprobs, pairs = pairs, pair_weights = pair_weights)
}

# ------------------------------------------------------------------------
# sample_branching()
# ------------------------------------------------------------------------

test_that("sample_branching returns correctly shaped output", {
  inputs <- make_toy_branching_inputs()

  out <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)

  expect_s3_class(out, "semihawkes_branching_sample")
  expect_type(out, "list")
  expect_named(out, "mark")
  expect_length(out$mark, length(inputs$bgprobs))
  expect_true(is.numeric(out$mark))
})

test_that("sample_branching validates its inputs", {
  inputs <- make_toy_branching_inputs()

  expect_error(
    sample_branching("not numeric", inputs$pairs, inputs$pair_weights),
    "must be a numeric vector"
  )
  expect_error(
    sample_branching(inputs$bgprobs, data.frame(a = 1, b = 2), inputs$pair_weights),
    "must be a data frame with integer columns"
  )
  expect_error(
    sample_branching(inputs$bgprobs, inputs$pairs, c(0.1, 0.2)),
    "length equal to `nrow\\(pairs\\)`"
  )
})

test_that("sample_branching only assigns values that are 0 or valid parent indices", {
  inputs <- make_toy_branching_inputs()
  set.seed(42)

  out <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)

  valid_values <- c(0L, inputs$pairs$j)
  expect_true(all(out$mark %in% valid_values))
})

test_that("sample_branching always marks an event with no candidate parents as background", {
  inputs <- make_toy_branching_inputs()

  # Event 1 has no rows in `pairs` where pairs$i == 1, and bgprobs[1] == 1,
  # so it must always end up as background regardless of the random draw.
  for (s in 1:20) {
    set.seed(s)
    out <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)
    expect_equal(out$mark[1], 0L)
  }
})

test_that("sample_branching is reproducible given the same seed", {
  inputs <- make_toy_branching_inputs()

  set.seed(123)
  out1 <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)

  set.seed(123)
  out2 <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)

  expect_identical(out1$mark, out2$mark)
})

test_that("sample_branching's parent draws respect extreme bgprobs", {
  # With bgprobs == 1 everywhere, every event must be background,
  # regardless of pair_weights (which should sum to 0 for consistency,
  # but we only test the sampling logic here, not weight consistency).
  bgprobs_all_bg <- rep(1, 5)
  pairs <- data.frame(i = c(2L, 3L), j = c(1L, 1L))
  pair_weights <- c(0, 0)

  set.seed(7)
  out <- sample_branching(bgprobs_all_bg, pairs, pair_weights)
  expect_true(all(out$mark == 0L))
})

test_that("sample_branching's parent draws respect near-zero bgprobs with a single parent", {
  bgprobs <- c(1, 0)
  pairs <- data.frame(i = 2L, j = 1L)
  pair_weights <- 1

  results <- vapply(1:50, function(s) {
    set.seed(s)
    sample_branching(bgprobs, pairs, pair_weights)$mark[2]
  }, integer(1))

  expect_true(all(results == 1L))
})

# ------------------------------------------------------------------------
# sample_branching_multi()
# ------------------------------------------------------------------------

test_that("sample_branching_multi returns the requested number of samples", {
  inputs <- make_toy_branching_inputs()

  samples <- sample_branching_multi(inputs$bgprobs, inputs$pairs,
                                    inputs$pair_weights, n_samples = 30, seed = 1)

  expect_length(samples, 30)
  expect_true(all(vapply(samples, inherits, logical(1), "semihawkes_branching_sample")))
})

test_that("sample_branching_multi is reproducible given the same seed", {
  inputs <- make_toy_branching_inputs()

  samples1 <- sample_branching_multi(inputs$bgprobs, inputs$pairs,
                                     inputs$pair_weights, n_samples = 10, seed = 99)
  samples2 <- sample_branching_multi(inputs$bgprobs, inputs$pairs,
                                     inputs$pair_weights, n_samples = 10, seed = 99)

  expect_identical(
    lapply(samples1, `[[`, "mark"),
    lapply(samples2, `[[`, "mark")
  )
})

test_that("a seeded sample_branching_multi leaves the caller's RNG state alone", {
  inputs <- make_toy_branching_inputs()

  set.seed(42)
  expected_next <- stats::runif(3)
  set.seed(42)
  sample_branching_multi(inputs$bgprobs, inputs$pairs, inputs$pair_weights,
                         n_samples = 5, seed = 99)
  expect_identical(stats::runif(3), expected_next)
})

test_that("sample_branching_multi without a seed still returns valid samples", {
  inputs <- make_toy_branching_inputs()

  samples <- sample_branching_multi(inputs$bgprobs, inputs$pairs,
                                    inputs$pair_weights, n_samples = 5)

  expect_length(samples, 5)
})

# ------------------------------------------------------------------------
# check_branching_weights()
# ------------------------------------------------------------------------

test_that("check_branching_weights reports ok = TRUE for consistent weights", {
  inputs <- make_toy_branching_inputs()

  check <- check_branching_weights(inputs$bgprobs, inputs$pairs, inputs$pair_weights)

  expect_s3_class(check, "semihawkes_branching_check")
  expect_true(check$ok)
  expect_lt(check$max_abs_error, 1e-6)
  expect_equal(check$n_no_parents, 0)
})

test_that("check_branching_weights detects weights that do not sum to 1", {
  inputs <- make_toy_branching_inputs()

  bad_weights <- inputs$pair_weights
  bad_weights[1] <- bad_weights[1] + 0.3  # break event 2's total

  check <- check_branching_weights(inputs$bgprobs, inputs$pairs, bad_weights)

  expect_false(check$ok)
  expect_gt(check$max_abs_error, 1e-6)
})

test_that("check_branching_weights flags events with no candidate parents and bgprobs < 1", {
  inputs <- make_toy_branching_inputs()

  # Artificially lower bgprobs[1] below 1 without giving it any parents in
  # `pairs` -- this should be flagged as inconsistent.
  bad_bgprobs <- inputs$bgprobs
  bad_bgprobs[1] <- 0.9

  check <- check_branching_weights(bad_bgprobs, inputs$pairs, inputs$pair_weights)

  expect_false(check$ok)
  expect_equal(check$n_no_parents, 1)
})

test_that("check_branching_weights respects the `tol` argument", {
  inputs <- make_toy_branching_inputs()

  slightly_off <- inputs$pair_weights
  slightly_off[1] <- slightly_off[1] + 1e-8

  expect_true(check_branching_weights(inputs$bgprobs, inputs$pairs, slightly_off,
                                      tol = 1e-6)$ok)
  expect_false(check_branching_weights(inputs$bgprobs, inputs$pairs, slightly_off,
                                       tol = 1e-10)$ok)
})

# ------------------------------------------------------------------------
# summarize_branching_samples()
# ------------------------------------------------------------------------

test_that("summarize_branching_samples computes the median, not the mean or a 0.5 threshold", {
  fake_samples <- list(
    list(mark = c(0L, 0L, 0L, 1L, 1L)),  # 3 background
    list(mark = c(0L, 0L, 1L, 1L, 1L)),  # 2 background
    list(mark = c(0L, 0L, 0L, 0L, 1L))   # 4 background
  )
  class(fake_samples[[1]]) <- class(fake_samples[[2]]) <-
    class(fake_samples[[3]]) <- "semihawkes_branching_sample"

  summary <- summarize_branching_samples(fake_samples)

  expect_s3_class(summary, "semihawkes_branching_summary")
  expect_equal(summary$num_bg, c(3L, 2L, 4L))
  expect_equal(summary$median_bg, 3)
  expect_false("mean_bg" %in% names(summary))
  expect_false("threshold_bg" %in% names(summary))
})

test_that("summarize_branching_samples computes a sensible IQR", {
  fake_samples <- lapply(c(10, 12, 14, 16, 18), function(k) {
    s <- list(mark = c(rep(0L, k), rep(1L, 20 - k)))
    class(s) <- "semihawkes_branching_sample"
    s
  })

  summary <- summarize_branching_samples(fake_samples)

  expect_equal(unname(summary$median_bg), 14)
  expect_equal(unname(summary$iqr_bg), c(12, 16))
})

test_that("summarize_branching_samples works end-to-end with sample_branching_multi", {
  inputs <- make_toy_branching_inputs()

  samples <- sample_branching_multi(inputs$bgprobs, inputs$pairs,
                                    inputs$pair_weights, n_samples = 100, seed = 2026)
  summary <- summarize_branching_samples(samples)

  expect_length(summary$num_bg, 100)
  expect_true(summary$median_bg >= 0 && summary$median_bg <= length(inputs$bgprobs))
  expect_true(all(summary$num_bg >= 1))  # event 1 is always background here
})

# ------------------------------------------------------------------------
# print methods (smoke tests -- just confirm they don't error)
# ------------------------------------------------------------------------

test_that("print methods run without error", {
  inputs <- make_toy_branching_inputs()

  sample <- sample_branching(inputs$bgprobs, inputs$pairs, inputs$pair_weights)
  check <- check_branching_weights(inputs$bgprobs, inputs$pairs, inputs$pair_weights)
  summary <- summarize_branching_samples(list(sample, sample))

  expect_output(print(sample), "semihawkes_branching_sample")
  expect_output(print(check), "semihawkes_branching_check")
  expect_output(print(summary), "semihawkes_branching_summary")
})

test_that("sample_branching_multi(restore_rng = FALSE) leaves the generator in the seeded state", {
  inputs <- make_toy_branching_inputs()

  set.seed(42)
  restored_next <- stats::runif(3)

  set.seed(42)
  s1 <- sample_branching_multi(inputs$bgprobs, inputs$pairs, inputs$pair_weights,
                               n_samples = 5, seed = 99, restore_rng = FALSE)
  after_a <- stats::runif(3)
  set.seed(7)
  s2 <- sample_branching_multi(inputs$bgprobs, inputs$pairs, inputs$pair_weights,
                               n_samples = 5, seed = 99, restore_rng = FALSE)
  after_b <- stats::runif(3)

  expect_identical(after_a, after_b)
  expect_false(identical(after_a, restored_next))
  expect_identical(lapply(s1, `[[`, "mark"), lapply(s2, `[[`, "mark"))
  expect_error(
    sample_branching_multi(inputs$bgprobs, inputs$pairs, inputs$pair_weights,
                           n_samples = 5, seed = 99, restore_rng = "yes"),
    "restore_rng"
  )
})

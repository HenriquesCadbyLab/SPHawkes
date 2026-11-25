# Ground-truth tests: does the estimator recover parameters it is given, and
# do the quantities it produces obey the invariants they must obey?
#
# See helper-hawkes-sim.R for why these exist alongside the reproduction
# tests in the rest of the suite.

# --- shared small fixture -----------------------------------------------
# Small and fast: enough events for the EM to run without degenerate
# histogram bins, small enough to fit in seconds rather than minutes.
fit_small <- function(seed = 42, n_iter = 4) {
  set.seed(seed)
  sim <- simulate_hawkes_branching(
    A = 0.4, n_background = 120, side = 8000, n_days = 400, sigma = 300, tau = 5
  )
  fit <- fit_semihawkes(
    sim$points, sim$region, date_col = "date",
    cell_size = 800, cache_dir = tempfile("spHawkes_recovery_"),
    bandwidth_args = list(k = 5, min_dist = 150, min_bandwidth = 300),
    n_iter = n_iter, verbose = FALSE
  )
  list(sim = sim, fit = fit)
}


test_that("the simulator produces a subcritical branching process with the requested shape", {
  set.seed(7)
  sim <- simulate_hawkes_branching(A = 0.4, n_background = 500)

  expect_s3_class(sim$points, "sf")
  expect_s3_class(sim$region, "sf")
  expect_identical(sim$truth$n_background, 500L)
  expect_gt(sim$truth$n_triggered, 0)

  # Total events for a subcritical process with immigration is n_bg / (1 - A)
  # in expectation, less whatever offspring fall outside the region or period.
  expect_lt(sim$truth$n, 500 / (1 - 0.4) * 1.2)
  expect_gt(sim$truth$n, 500)

  # Every event must lie inside the region and the study period.
  expect_true(all(sim$events$day >= 0 & sim$events$day < 730))
  inside <- sf::st_intersects(sim$points, sim$region, sparse = FALSE)[, 1]
  expect_true(all(inside))

  # Parents strictly precede their offspring (no same-day pairs by default).
  expect_true(all(sim$events$generation >= 0L))
})


test_that("fitted quantities obey their invariants at every iteration", {
  out <- fit_small()
  fit <- out$fit

  # Background probabilities are probabilities.
  expect_true(all(fit$bgprobs >= 0 & fit$bgprobs <= 1),
              info = "bgprobs outside [0, 1]")
  expect_false(anyNA(fit$bgprobs))

  # Scalar parameters are finite and positive.
  expect_true(is.finite(fit$A) && fit$A > 0)
  expect_true(is.finite(fit$mu) && fit$mu > 0)
  expect_true(is.finite(fit$llik))

  # The intensity is strictly positive wherever an event was observed --
  # otherwise the log-likelihood would be undefined.
  expect_true(all(is.finite(fit$lambda_at_events)))
  expect_true(all(fit$lambda_at_events > 0))

  # history is one row per iteration, starting at Loop 0, and its last row
  # agrees with the object's own final scalars.
  expect_equal(fit$history$iteration, seq(0, fit$n_iter))
  last <- fit$history[nrow(fit$history), ]
  expect_equal(last$mu, fit$mu)
  expect_equal(last$A, fit$A)
  expect_equal(last$llik, fit$llik)
})


test_that("re-estimated components keep their normalisations", {
  out <- fit_small()
  fit <- out$fit

  # The background multipliers are identified only up to scale, and are
  # normalised to mean 1 so that mu carries the level.
  expect_equal(mean(fit$weekly$basevalue), 1, tolerance = 1e-6)
  expect_equal(mean(fit$monthly$basevalue), 1, tolerance = 1e-6)
  expect_equal(mean(fit$trend$basevalue), 1, tolerance = 1e-6)

  # The spatial background is mean-normalised over in-region cells.
  in_region <- fit$background$marks > 0
  expect_equal(mean(fit$background$basevalue[in_region]), 1, tolerance = 1e-6)

  # Both triggering kernels are densities: they integrate to 1, so that A
  # alone carries the expected number of offspring.
  g_step <- fit$temporal_trigger$base[2] - fit$temporal_trigger$base[1]
  expect_equal(simpson(fit$temporal_trigger$basevalue, g_step), 1, tolerance = 1e-6)

  h_dx <- fit$spatial_trigger$x[2] - fit$spatial_trigger$x[1]
  h_dy <- fit$spatial_trigger$y[2] - fit$spatial_trigger$y[1]
  expect_equal(simpson.2D(fit$spatial_trigger$basevalue, dx = h_dx, dy = h_dy), 1,
               tolerance = 1e-6)
})


test_that("the EM log-likelihood never decreases across iterations", {
  out <- fit_small(n_iter = 6)
  deltas <- diff(out$fit$history$llik)

  # This is the defining property of EM: each iteration must not decrease the
  # observed-data log-likelihood. A negative step means the M-step is not
  # maximising what the E-step assumed -- the single most informative check
  # that the update formulas are mutually consistent.
  expect_true(all(deltas >= 0),
              info = paste0("log-likelihood decreased at iteration(s) ",
                            paste(which(deltas < 0), collapse = ", "),
                            "; deltas: ", paste(round(deltas, 6), collapse = ", ")))

  # And it should be converging, not wandering: later steps smaller than the first.
  expect_lt(abs(deltas[length(deltas)]), abs(deltas[1]))
})


test_that("the fitted intensity integrates to exactly the number of events", {
  out <- fit_small()
  n_events <- nrow(out$fit$events$events)

  # The M-step sets mu so that mu * bg_at_all + A * trig_at_all == n exactly.
  # It is an algebraic identity rather than an approximation, so it holds to
  # floating-point precision and is a sharp check on that update.
  expect_equal(out$fit$lambda_at_all, n_events, tolerance = 1e-8)
  expect_true(all(abs(out$fit$history$lambda_at_all - n_events) < 1e-6))
})


# --- parameter recovery from a known process ----------------------------
#
# This is the test the rest of the suite cannot provide: it generates data
# from a process whose parameters are known, fits the model, and asks whether
# the estimates come back.
#
# It is opt-in (see `skip_if_not_slow()`) because it fits a full model, and it
# runs on one CI job rather than all four.
#
# MEASURED BEHAVIOUR, spHawkes 0.1.1. Eight fits of
# `simulate_hawkes_branching(A = 0.4)` with its other defaults (seeds 1-8,
# n = 629-687 events), `cell_size = 750`, 30 iterations, `tol = 1e-4`:
#
#                    4 km reach                default 1 km reach
#   A_hat            mean 0.61, +54%           mean 0.47, +17%
#                    (+40% to +76% per fit)    (+8% to +24% per fit)
#   triggered share  +13.3 points too high     +2.3 points too high
#   converged        0 of 8 (all max_iter)     0 of 8 (all max_iter)
#
# Most of the excess comes from the triggering kernel's reach. Offspring are
# displaced from their parent by a Gaussian with sd 300 m per coordinate, so
# only 0.4% of them land more than 1 km away. At either reach, the triggering
# weight on pairs of events less than 1 km apart adds up to within 5% of the
# true number of triggered events. With a 4 km reach, though, pairs more than
# 1 km apart add another 80 events' worth of weight per fit on average (63 to
# 108): the kernel's tail picks up pairs of unrelated events. The default
# 1 km reach cuts that to 4. The excess that remains at 1 km, larger in A
# than in the share, is not yet attributed.
#
# Also measured: on a background that is homogeneous in time, the fitted
# trend, monthly and weekly multipliers should come out flat at 1, but their
# values span about 0.5, 0.24 and 0.15 on average, at either reach. They
# absorb noise as structure.
#
# This test fits with the default 1 km reach. Its bounds are deliberately
# wide: they are set to catch the estimator breaking, not to certify it as
# unbiased.

test_that("A and the triggering kernels are recovered from a simulated process", {
  # NOTE: deliberately no skip_on_cran(). Under `R CMD check` the NOT_CRAN
  # variable is unset, so skip_on_cran() would skip this test on CI too --
  # exactly where it is meant to run. The SPHAWKES_SLOW_TESTS gate above
  # already keeps it off CRAN, since CRAN will never set it.
  skip_if_not_slow()

  A_true <- 0.4
  set.seed(2024)
  sim <- simulate_hawkes_branching(
    A = A_true, n_background = 250, side = 12000, n_days = 550,
    sigma = 300, tau = 5
  )

  fit <- fit_semihawkes(
    sim$points, sim$region, date_col = "date",
    cell_size = 750, cache_dir = tempfile("spHawkes_recovery_slow_"),
    bandwidth_args = list(k = 5, min_dist = 150, min_bandwidth = 300),
    n_iter = 20, verbose = FALSE
  )

  # --- the branching ratio ---------------------------------------------
  # Wide by design; see the note above. Measured 0.46 for this seed.
  expect_gt(fit$A, 0.5 * A_true)
  expect_lt(fit$A, 2.2 * A_true)

  # --- the share of events attributed to triggering ---------------------
  # Recovered much better than A itself. Measured |diff| = 0.009 for this
  # seed.
  share_hat <- 1 - mean(fit$bgprobs)
  expect_lt(abs(share_hat - sim$truth$triggered_share), 0.15,
            label = sprintf("|share_hat - share_true| (hat %.3f, true %.3f)",
                            share_hat, sim$truth$triggered_share))

  # --- shape of the temporal triggering kernel --------------------------
  # True lags are Exp(mean 5 days). The estimate is smoothed with a 5-day
  # bandwidth so it is wider than the truth, but its mass must still sit at
  # short lags rather than spread over the 90-day window. Measured 0.93.
  g <- fit$temporal_trigger
  g_step <- g$base[2] - g$base[1]
  mass_20d <- simpson(ifelse(g$base <= 20, g$basevalue, 0), g_step)
  expect_gt(mass_20d, 0.75)

  # --- shape of the spatial triggering kernel ---------------------------
  # True offspring displacement is Gaussian with sd 300 m per coordinate.
  # The estimate is smoothed at 250 m on a 200 m grid, so it is markedly more
  # diffuse than the truth, but must remain concentrated near the origin of
  # the +/-1 km grid it lives on. Measured 0.56 within 500 m, against 0.75
  # for the true kernel and 0.23 for one spread evenly over the grid's disc.
  h <- fit$spatial_trigger
  hx <- h$x[2] - h$x[1]
  hy <- h$y[2] - h$y[1]
  radius <- sqrt((h$x %o% rep(1, length(h$y)))^2 +
                   (rep(1, length(h$x)) %o% h$y)^2)
  mass_500m <- simpson.2D(h$basevalue * (radius <= 500), dx = hx, dy = hy)
  expect_gt(mass_500m, 0.40)

  # The kernel must peak at the origin, not somewhere else.
  peak <- which(h$basevalue == max(h$basevalue), arr.ind = TRUE)[1, ]
  expect_lt(abs(h$x[peak[["row"]]]), 750)
  expect_lt(abs(h$y[peak[["col"]]]), 750)
})

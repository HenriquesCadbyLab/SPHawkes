# Simulate a spatio-temporal Hawkes process with KNOWN parameters.
#
# This exists so the estimator can be checked against a ground truth rather
# than against itself. The rest of the suite is largely reproduction testing:
# it re-runs the same exported building blocks a second way and compares, so
# it confirms that `em_step()` composes its parts correctly, but would pass
# just as happily if a shared formula were wrong. This helper closes that gap.
#
# The generating process mirrors the model `fit_semihawkes()` estimates:
#
#   background  homogeneous Poisson over the region and study period, so the
#               trend / weekly / monthly components have no real signal to
#               find and should come out flat at about 1
#   offspring   every event independently spawns Poisson(A) children, each
#               Exp(tau) days later and Gaussian(sigma) away in each
#               coordinate, cascading to extinction
#
# Because g and h are normalised to integrate to 1, `A` is exactly the
# expected number of direct offspring per event -- the quantity
# `fit_semihawkes()` reports as `$A`. For a subcritical process it is also
# the long-run share of events that are offspring, which is what
# `1 - mean(fit$bgprobs)` estimates.
#
# NOTE ON TIES: lags are forced to at least one whole day. The package
# evaluates its temporal triggering kernel at `g(0) = 0` exactly, so
# same-day pairs contribute nothing and would bias `A` downwards for reasons
# unrelated to whatever is being tested. Keeping offspring at least a day
# after their parent isolates that out. Set `min_lag_days = 0` to include
# same-day offspring and see the effect.

simulate_hawkes_branching <- function(A = 0.4,
                                      n_background = 400,
                                      side = 15000,
                                      n_days = 730,
                                      sigma = 300,
                                      tau = 5,
                                      start_date = as.Date("2021-01-03"),
                                      crs = 32717,
                                      min_lag_days = 1,
                                      max_generations = 12L) {
  stopifnot(A >= 0, n_background >= 1, side > 0, n_days > 1, sigma > 0, tau > 0)

  events <- data.frame(
    x   = stats::runif(n_background, 0, side),
    y   = stats::runif(n_background, 0, side),
    day = sample.int(n_days, n_background, replace = TRUE) - 1L,
    generation = 0L
  )

  parents <- events
  generation <- 0L
  while (nrow(parents) > 0 && generation < max_generations) {
    generation <- generation + 1L
    n_kids <- stats::rpois(nrow(parents), A)
    if (sum(n_kids) == 0) break

    idx <- rep(seq_len(nrow(parents)), n_kids)
    lag <- stats::rexp(length(idx), rate = 1 / tau)
    lag <- pmax(min_lag_days, ceiling(lag))

    kids <- data.frame(
      x   = parents$x[idx] + stats::rnorm(length(idx), 0, sigma),
      y   = parents$y[idx] + stats::rnorm(length(idx), 0, sigma),
      day = parents$day[idx] + lag,
      generation = generation
    )
    # Offspring leaving the region or the study period are unobserved, exactly
    # as they would be in real data.
    kids <- kids[kids$x > 0 & kids$x < side &
                   kids$y > 0 & kids$y < side &
                   kids$day < n_days, , drop = FALSE]
    if (nrow(kids) == 0) break

    events <- rbind(events, kids)
    parents <- kids
  }

  events <- events[order(events$day), ]
  rownames(events) <- NULL

  square <- sf::st_polygon(list(matrix(
    c(0, 0, side, 0, side, side, 0, side, 0, 0), ncol = 2, byrow = TRUE
  )))

  list(
    points = sf::st_as_sf(
      data.frame(date = start_date + events$day, x = events$x, y = events$y),
      coords = c("x", "y"), crs = crs
    ),
    region = sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = crs)),
    events = events,
    truth  = list(
      A = A, sigma = sigma, tau = tau,
      n = nrow(events),
      n_background = sum(events$generation == 0L),
      n_triggered  = sum(events$generation > 0L),
      triggered_share = mean(events$generation > 0L)
    )
  )
}

# The full recovery fit takes a few minutes, which is too slow to run on
# every push across four platforms. It is therefore opt-in: set
#   SPHAWKES_SLOW_TESTS=true
# to run it locally, and one CI job sets it so the check still happens on
# every pull request without multiplying the cost by four.
skip_if_not_slow <- function() {
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("SPHAWKES_SLOW_TESTS")), "true"),
    "slow recovery test; set SPHAWKES_SLOW_TESTS=true to run"
  )
}

# Generates the `sim_events` / `sim_region` example datasets shipped in
# data/. Run with `source("data-raw/make-example-data.R")` from the package
# root to regenerate them.
#
# The pattern is SIMULATED, not real crime data. It is drawn from a
# branching process deliberately shaped like the model this package fits, so
# every fitted component has something real to recover:
#
#   * spatial background : three Gaussian hotspots plus a uniform floor
#   * long-term trend    : immigration rate rising over the three years
#   * day of week        : weekends elevated relative to weekdays
#   * triggering         : each event spawns Poisson(A) offspring, close in
#                          space (~300 m) and soon after in time (~5 days)
#
# Known generating values, useful as a target for a recovery test:
#   A (branching ratio)  = 0.35
#   trigger length scale = 300 m, trigger time scale = 5 days

library(sf)

set.seed(20260818)

## ---- CRS ---------------------------------------------------------------
# UTM zone 17S (metres) -- the same projection family as the Guayaquil
# application, so the package's metre-based defaults are sensible.
#
# EPSG:32717's full WKT carries degree signs in its USAGE[AREA[...]]
# description, which trips `R CMD check`'s "data for non-ASCII characters"
# warning when an sf object using it is saved to data/. That block is
# descriptive metadata only, so it is stripped here: the result still
# compares equal to `st_crs(32717)`, keeps the EPSG code, and reports
# `units = "m"`.
utm17s <- local({
  wkt <- sf::st_crs(32717)$wkt
  crs <- sf::st_crs(sub("\\s*USAGE\\[.*?\\]\\],\\s*", "\n    ", wkt))
  stopifnot(crs == sf::st_crs(32717), identical(crs$units, "m"))
  crs
})

## ---- study region -----------------------------------------------------
side <- 8000
region_poly <- st_polygon(list(matrix(
  c(0, 0,
    side, 0,
    side, side,
    0, side,
    0, 0),
  ncol = 2, byrow = TRUE
)))
sim_region <- st_sf(
  name = "Example city",
  geometry = st_sfc(region_poly, crs = utm17s)
)

## ---- background (immigrant) events ------------------------------------
start_date <- as.Date("2021-01-01")
end_date   <- as.Date("2024-01-01")
n_days     <- as.integer(end_date - start_date)

hotspots <- data.frame(
  x = c(2000, 6000, 4200),
  y = c(6200, 2600, 4000),
  sd = c(700, 900, 500),
  w = c(0.34, 0.30, 0.18)
)
p_uniform <- 1 - sum(hotspots$w)

n_background <- 260

# Time: rising trend x weekend effect, sampled by inversion over days.
day_index <- seq_len(n_days) - 1L
trend    <- 0.55 + 0.9 * (day_index / n_days)^1.6
wday     <- as.POSIXlt(start_date + day_index)$wday
weekend  <- ifelse(wday %in% c(0, 6), 1.55, 0.85)
day_prob <- trend * weekend
day_prob <- day_prob / sum(day_prob)

bg_days <- sample(day_index, n_background, replace = TRUE, prob = day_prob)

# Space: mixture of hotspots and a uniform floor, rejection-sampled to the
# region so nothing lands outside it.
draw_bg_xy <- function(n) {
  out <- matrix(NA_real_, nrow = 0, ncol = 2)
  while (nrow(out) < n) {
    m <- n - nrow(out)
    comp <- sample(
      c(seq_len(nrow(hotspots)), 0L), m, replace = TRUE,
      prob = c(hotspots$w, p_uniform)
    )
    xy <- matrix(NA_real_, nrow = m, ncol = 2)
    is_unif <- comp == 0L
    if (any(is_unif)) {
      xy[is_unif, ] <- cbind(runif(sum(is_unif), 0, side),
                             runif(sum(is_unif), 0, side))
    }
    for (k in seq_len(nrow(hotspots))) {
      sel <- comp == k
      if (any(sel)) {
        xy[sel, ] <- cbind(rnorm(sum(sel), hotspots$x[k], hotspots$sd[k]),
                           rnorm(sum(sel), hotspots$y[k], hotspots$sd[k]))
      }
    }
    keep <- xy[, 1] > 0 & xy[, 1] < side & xy[, 2] > 0 & xy[, 2] < side
    out <- rbind(out, xy[keep, , drop = FALSE])
  }
  out[seq_len(n), , drop = FALSE]
}

bg_xy <- draw_bg_xy(n_background)

events <- data.frame(
  x = bg_xy[, 1],
  y = bg_xy[, 2],
  day = bg_days,
  generation = 0L
)

## ---- offspring (triggered) events -------------------------------------
A            <- 0.35   # expected direct offspring per event
space_scale  <- 300    # metres
time_scale   <- 5      # days

parents <- events
generation <- 0L
while (nrow(parents) > 0 && generation < 8L) {
  generation <- generation + 1L
  n_kids <- rpois(nrow(parents), A)
  if (sum(n_kids) == 0) break

  idx <- rep(seq_len(nrow(parents)), n_kids)
  kids <- data.frame(
    x = parents$x[idx] + rnorm(length(idx), 0, space_scale),
    y = parents$y[idx] + rnorm(length(idx), 0, space_scale),
    day = parents$day[idx] + ceiling(rexp(length(idx), rate = 1 / time_scale)),
    generation = generation
  )

  kids <- kids[kids$x > 0 & kids$x < side &
                 kids$y > 0 & kids$y < side &
                 kids$day < n_days, , drop = FALSE]
  if (nrow(kids) == 0) break

  events <- rbind(events, kids)
  parents <- kids
}

events <- events[order(events$day), ]
rownames(events) <- NULL

## ---- assemble the sf layer --------------------------------------------
sim_events <- st_as_sf(
  data.frame(
    id   = seq_len(nrow(events)),
    date = start_date + events$day,
    x    = events$x,
    y    = events$y
  ),
  coords = c("x", "y"),
  crs = utm17s
)

message(sprintf(
  "simulated %d events (%d background, %d triggered = %.1f%%) over %s to %s",
  nrow(events), sum(events$generation == 0L), sum(events$generation > 0L),
  100 * mean(events$generation > 0L), start_date, max(sim_events$date)
))

usethis::use_data(sim_events, sim_region, overwrite = TRUE)

# Small synthetic dataset shared by tests: a 1000m x 1000m square study
# region (projected CRS, EPSG:32717 just as a stand-in metric CRS) with
# uniformly random points and dates, plus a few duplicated/near-duplicate
# points to exercise the adaptive bandwidth fallback logic.

make_synthetic_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  square <- sf::st_polygon(list(matrix(
    c(0, 0, 1000, 0, 1000, 1000, 0, 1000, 0, 0),
    ncol = 2, byrow = TRUE
  )))
  region <- sf::st_sf(id = 1, geometry = sf::st_sfc(square, crs = 32717))

  x <- stats::runif(n, 0, 1000)
  y <- stats::runif(n, 0, 1000)
  dates <- as.Date("2020-01-01") + sample(0:364, n, replace = TRUE)

  pts <- sf::st_as_sf(
    data.frame(x = x, y = y, event_date = dates),
    coords = c("x", "y"), crs = 32717
  )

  list(points = pts, region = region)
}

# Lightweight fitted model used by the goodness-of-fit tests.
#
# Uses `n_iter = 0`, so it returns the initial Loop-0 fit without running any
# EM update. This is sufficient for testing the goodness-of-fit functions and
# avoids the extra cost and fragility of fitting synthetic data.
make_gof_fit <- function(n = 250, seed = 21) {
  synth <- make_synthetic_data(n = n, seed = seed)

  cell_size <- 150
  padded_region <- sf::st_buffer(
    synth$region,
    dist = 4 * cell_size
  )

  fit_semihawkes(
    points = synth$points,
    region = padded_region,
    date_col = "event_date",
    cell_size = cell_size,
    cache_dir = tempfile("spHawkes-gof-marks-"),
    bandwidth_args = list(
      k = 3,
      min_dist = 10,
      min_bandwidth = 50
    ),
    n_iter = 0,
    verbose = FALSE
  )
}

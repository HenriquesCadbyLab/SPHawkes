temp_mark_dir <- function() {
  d <- tempfile("marks_")
  dir.create(d)
  d
}

make_ready_events3 <- function(n = 30, seed = 11) {
  synth <- make_synthetic_data(n = n, seed = seed)
  ev <- prepare_events(synth$points, synth$region, date_col = "event_date",
                        start_date = "2020-01-01")
  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
  ev
}

test_that("trigger_pairs matches a brute-force double loop", {
  ev <- make_ready_events3(n = 25)
  pairs <- trigger_pairs(ev, window_days = 90, window_dist = 4000)

  d <- ev$events
  n <- nrow(d)
  brute <- data.frame(i = integer(0), j = integer(0))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (d$days[i] > d$days[j] &&
          d$days[i] <= d$days[j] + 90 &&
          abs(d$coorx[i] - d$coorx[j]) <= 4000 &&
          abs(d$coory[i] - d$coory[j]) <= 4000) {
        brute <- rbind(brute, data.frame(i = i, j = j))
      }
    }
  }

  expect_equal(nrow(pairs), nrow(brute))
  expect_true(nrow(pairs) > 0)
  # order-independent comparison
  key_a <- paste(pairs$i, pairs$j)
  key_b <- paste(brute$i, brute$j)
  expect_equal(sort(key_a), sort(key_b))
  # sanity: i == j never appears (strict day comparison already excludes it)
  expect_true(all(pairs$i != pairs$j))
})

test_that("trigger_pairs respects window_days and window_dist", {
  ev <- make_ready_events3(n = 25)
  wide <- trigger_pairs(ev, window_days = 90, window_dist = 4000)
  narrow <- trigger_pairs(ev, window_days = 10, window_dist = 4000)
  expect_true(nrow(narrow) <= nrow(wide))

  d <- ev$events
  if (nrow(narrow) > 0) {
    expect_true(all(d$days[narrow$i] - d$days[narrow$j] <= 10))
  }
})

test_that("trigger_pairs returns exactly what comparing every pair returns, in the same order", {
  # The all-pairs version trigger_pairs() used before, kept as the reference.
  all_pairs <- function(ev, window_days, window_dist) {
    ev <- ev$events
    day_after <- outer(ev$days, ev$days, `>`)
    day_window <- outer(ev$days, ev$days, function(di, dj) di <= dj + window_days)
    x_window <- outer(ev$coorx, ev$coorx, function(xi, xj) abs(xi - xj) <= window_dist)
    y_window <- outer(ev$coory, ev$coory, function(yi, yj) abs(yi - yj) <= window_dist)
    idx <- which(day_after & day_window & x_window & y_window, arr.ind = TRUE)
    data.frame(i = as.integer(idx[, "row"]), j = as.integer(idx[, "col"]))
  }

  ev <- make_ready_events3(n = 120, seed = 4) # one year: many events share a day
  expect_gt(anyDuplicated(ev$events$days), 0)
  for (w in list(c(90, 4000), c(10, 300), c(1, 50), c(0, 4000))) {
    expect_identical(trigger_pairs(ev, window_days = w[1], window_dist = w[2]),
                     all_pairs(ev, w[1], w[2]))
  }
  expect_identical(nrow(trigger_pairs(ev, window_days = 0, window_dist = 4000)), 0L)
})

test_that("trigger_spatial_marks caches marks matching a brute-force sf check", {
  ev <- make_ready_events3(n = 8)
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()

  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)

  # brute force for every event: shift the grid by the event's coords and
  # test inside-region membership directly via sf
  d <- ev$events
  nx <- length(ts$x); ny <- length(ts$y)
  X <- ts$x %o% rep(1, ny)
  Y <- rep(1, nx) %o% ts$y
  for (i in seq_len(nrow(d))) {
    pts <- sf::st_as_sf(data.frame(x = as.vector(X) + d$coorx[i], y = as.vector(Y) + d$coory[i]),
                        coords = c("x", "y"), crs = sf::st_crs(ev$region))
    brute_mark <- matrix(lengths(sf::st_intersects(pts, ev$region)) > 0, nrow = nx, ncol = ny)
    expect_identical(readRDS(file.path(cache_dir, sprintf("mark_%06d.rds", i))), brute_mark)
  }
  expect_true(all(file.exists(file.path(cache_dir, sprintf("mark_%06d.rds", seq_len(nrow(d)))))))
})

test_that("trigger_spatial_marks counts grid cells exactly on the boundary as inside, as sf does", {
  square <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32717))
  # Round coordinates and a 200 m kernel grid put many shifted grid cells
  # exactly on the square's edges and corners.
  pts <- sf::st_as_sf(data.frame(x = c(200, 400, 600, 800, 200, 600),
                                 y = c(200, 400, 800, 600, 600, 200),
                                 event_date = as.Date("2020-01-01") + 0:5),
                      coords = c("x", "y"), crs = 32717)
  ev <- prepare_events(pts, square, date_col = "event_date")
  ts <- trigger_spatial_init(x_range = 1000, x_step = 200)
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)

  nx <- length(ts$x); ny <- length(ts$y)
  X <- as.vector(ts$x %o% rep(1, ny)); Y <- as.vector(rep(1, nx) %o% ts$y)
  on_edge <- 0
  for (i in seq_len(nrow(ev$events))) {
    x <- X + ev$events$coorx[i]; y <- Y + ev$events$coory[i]
    on_edge <- on_edge + sum((x %in% c(0, 1000) & y >= 0 & y <= 1000) |
                               (y %in% c(0, 1000) & x >= 0 & x <= 1000))
    p <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = 32717)
    expect_identical(readRDS(file.path(cache_dir, sprintf("mark_%06d.rds", i))),
                     matrix(lengths(sf::st_intersects(p, ev$region)) > 0, nrow = nx, ncol = ny))
  }
  expect_gt(on_edge, 0)
})

test_that("the fast point-in-region test agrees with sf, including on and next to the boundary", {
  # An L shape, so the boundary has both outward and inward corners.
  L <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 400), c(400, 400), c(400, 1000), c(0, 1000), c(0, 0)))),
    crs = 32717))
  on_boundary <- rbind(c(0, 0), c(500, 0), c(1000, 200), c(1000, 400), c(700, 400),
                       c(400, 400), c(400, 700), c(400, 1000), c(200, 1000), c(0, 500))
  next_to_it <- rbind(on_boundary + 1e-7, on_boundary - 1e-7,
                      cbind(on_boundary[, 1] + 1e-7, on_boundary[, 2] - 1e-7))
  lattice <- as.matrix(expand.grid(seq(-50, 1050, 25), seq(-50, 1050, 25)))
  pts <- rbind(on_boundary, next_to_it, lattice)

  expected <- lengths(sf::st_intersects(
    sf::st_as_sf(data.frame(x = pts[, 1], y = pts[, 2]), coords = c("x", "y"), crs = 32717),
    L)) > 0
  in_region <- .region_point_test(L, eps = 0.2)
  expect_identical(in_region(pts[, 1], pts[, 2]), expected)
  expect_true(all(in_region(on_boundary[, 1], on_boundary[, 2])))

  # A region thinner than the band: every point is tested with sf instead.
  thin <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 0.1), c(0, 0.1), c(0, 0)))), crs = 32717))
  x <- c(500, 500, 500, -1); y <- c(0.05, 0.1, 0.2, 0.05)
  expect_identical(.region_point_test(thin, eps = 0.2)(x, y), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("trigger_spatial_marks skips existing files unless force = TRUE", {
  ev <- make_ready_events3(n = 5)
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()

  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
  fn <- file.path(cache_dir, sprintf("mark_%06d.rds", 1))
  mtime1 <- file.info(fn)$mtime

  Sys.sleep(1.1)
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir, force = FALSE)
  expect_equal(file.info(fn)$mtime, mtime1)

  trigger_spatial_marks(ev, ts, cache_dir = cache_dir, force = TRUE)
  expect_true(file.info(fn)$mtime > mtime1)
})

test_that("trigger_spatial_marks nprocs = 1 and nprocs = 2 produce identical caches", {
  ev <- make_ready_events3(n = 10)
  ts <- trigger_spatial_init()
  dir1 <- temp_mark_dir()
  dir2 <- temp_mark_dir()

  trigger_spatial_marks(ev, ts, cache_dir = dir1, nprocs = 1)
  trigger_spatial_marks(ev, ts, cache_dir = dir2, nprocs = 2)

  for (i in seq_len(nrow(ev$events))) {
    expect_equal(readRDS(file.path(dir1, sprintf("mark_%06d.rds", i))),
                 readRDS(file.path(dir2, sprintf("mark_%06d.rds", i))))
  }

  # A cluster passed in is used and left running.
  dir3 <- temp_mark_dir()
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  trigger_spatial_marks(ev, ts, cache_dir = dir3, nprocs = cl)
  for (i in seq_len(nrow(ev$events))) {
    expect_equal(readRDS(file.path(dir3, sprintf("mark_%06d.rds", i))),
                 readRDS(file.path(dir1, sprintf("mark_%06d.rds", i))))
  }
  expect_equal(unlist(parallel::clusterEvalQ(cl, 1 + 1)), c(2, 2))
})

test_that("trigger_edge_correction matches a brute-force per-event computation", {
  ev <- make_ready_events3(n = 12)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)

  ec <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)
  expect_s3_class(ec, "semihawkes_trigger_edge_correction")
  expect_length(ec$temporal_edge_correction, nrow(ev$events))
  expect_length(ec$spatial_edge_correction, nrow(ev$events))

  d <- ev$events
  step <- tt$base[2] - tt$base[1]
  brute_temporal_1 <- sum(tt$fun(seq(0, ev$TT - d$days[1], step) + 6e-6))
  expect_equal(ec$temporal_edge_correction[1], brute_temporal_1, tolerance = 1e-8)

  mark1 <- readRDS(file.path(cache_dir, sprintf("mark_%06d.rds", 1)))
  dx <- ts$x[2] - ts$x[1]; dy <- ts$y[2] - ts$y[1]
  brute_spatial_1 <- simpson.2D(mark1 * ts$basevalue, dx = dx, dy = dy)
  expect_equal(ec$spatial_edge_correction[1], brute_spatial_1, tolerance = 1e-8)

  # temporal_repetance_fun should equal exact counts at the grid nodes
  brute_repetance <- vapply(tt$base, function(b) sum(b < ev$TT - d$days), numeric(1))
  expect_equal(ec$temporal_repetance_fun(tt$base), brute_repetance, tolerance = 1e-8)
})

test_that("trigger_edge_correction's temporal factors equal sums over each event's whole remaining window", {
  ev <- make_ready_events3(n = 20)
  tt <- trigger_temporal_init() # reaches 90 days, in a one-year window
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)

  ec <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)

  step <- tt$base[2] - tt$base[1]
  brute <- vapply(ev$events$days, function(d) {
    sum(tt$fun(seq(0, ev$TT - d, step) + 6e-6))
  }, numeric(1))
  expect_identical(ec$temporal_edge_correction, brute)
})

test_that("trigger_edge_correction gives the same result with the marks read in advance", {
  ev <- make_ready_events3(n = 15)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)

  from_disk <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)
  marks <- .read_spatial_marks(cache_dir, nrow(ev$events))
  in_memory <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir, marks = marks)

  for (nm in c("temporal_edge_correction", "spatial_edge_correction",
               "temporal_repetance", "spatial_repetance")) {
    expect_identical(in_memory[[nm]], from_disk[[nm]], label = nm)
  }
  expect_error(trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir,
                                       marks = marks[-1]),
               "`marks` does not match")
})

test_that("trigger_edge_correction errors when marks are missing", {
  ev <- make_ready_events3(n = 5)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  empty_dir <- temp_mark_dir()
  expect_error(trigger_edge_correction(ev, tt, ts, cache_dir = empty_dir),
               "trigger_spatial_marks")
})

test_that("trigger_pair_weights matches the manual formula", {
  ev <- make_ready_events3(n = 15)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  pairs <- trigger_pairs(ev)
  d <- ev$events
  lambda_at_events <- rep(1, nrow(d)) # dummy for the formula check
  A <- 0.3

  w <- trigger_pair_weights(ev, pairs, A, lambda_at_events, tt, ts)

  brute <- A * tt$fun(d$days[pairs$i] - d$days[pairs$j]) *
    ts$fun(d$coorx[pairs$i] - d$coorx[pairs$j], d$coory[pairs$i] - d$coory[pairs$j]) /
    lambda_at_events[pairs$i]

  expect_equal(w, brute, tolerance = 1e-10)
})

test_that("trigger_reestimate_temporal returns a normalized, drop-in trigger object", {
  ev <- make_ready_events3(n = 25)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
  ec <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)

  pairs <- trigger_pairs(ev)
  skip_if(nrow(pairs) == 0, "no triggering pairs in this synthetic draw")

  lambda_at_events <- rep(1, nrow(ev$events))
  w <- trigger_pair_weights(ev, pairs, A = 0.3, lambda_at_events, tt, ts)

  tt2 <- trigger_reestimate_temporal(ev, pairs, w, ec, tt)
  expect_s3_class(tt2, "semihawkes_trigger_temporal")
  expect_equal(tt2$base, tt$base)
  step <- tt$base[2] - tt$base[1]
  expect_equal(simpson(tt2$basevalue, step), 1, tolerance = 1e-6)
  expect_equal(tt2$fun(-5), 0)
  expect_equal(tt2$fun(1000), 0)
})

test_that("trigger_reestimate_spatial returns a normalized, drop-in trigger object", {
  ev <- make_ready_events3(n = 25)
  tt <- trigger_temporal_init()
  # A grid sized to the 1000x1000 synthetic region, so that some pair can
  # realise most offsets on it (a much wider grid drives the spatial
  # repetance surface to zero almost everywhere). This keeps the test in the
  # well-behaved regime; the asymmetric-grid test below deliberately probes
  # the y-mesh fix instead of normalization.
  ts <- trigger_spatial_init(x_range = 500, y_range = 500, x_step = 50, y_step = 50)
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
  ec <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)

  # window_dist matched to the grid's range: pair offsets outside the
  # trigger grid are genuinely out-of-range for hist.weighted.2D() (which
  # deliberately does not silently coerce them, see its docs), and this
  # small region's real pair offsets otherwise exceed a +-500 grid.
  pairs <- trigger_pairs(ev, window_dist = 500)
  skip_if(nrow(pairs) == 0, "no triggering pairs in this synthetic draw")

  lambda_at_events <- rep(1, nrow(ev$events))
  w <- trigger_pair_weights(ev, pairs, A = 0.3, lambda_at_events, tt, ts)

  ts2 <- trigger_reestimate_spatial(ev, pairs, w, ec, ts)
  expect_s3_class(ts2, "semihawkes_trigger_spatial")
  expect_equal(ts2$x, ts$x)
  expect_equal(ts2$y, ts$y)
  dx <- ts$x[2] - ts$x[1]; dy <- ts$y[2] - ts$y[1]
  expect_equal(simpson.2D(ts2$basevalue, dx = dx, dy = dy), 1, tolerance = 1e-6)
  expect_equal(ts2$fun(1e7, 1e7), 0)
})

test_that("trigger_reestimate_spatial's y-mesh fix actually changes the result on an asymmetric grid", {
  ev <- make_ready_events3(n = 25)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init(x_range = 4000, y_range = 4000, x_step = 200, y_step = 400)
  cache_dir <- temp_mark_dir()
  trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
  ec <- trigger_edge_correction(ev, tt, ts, cache_dir = cache_dir)

  pairs <- trigger_pairs(ev)
  skip_if(nrow(pairs) == 0, "no triggering pairs in this synthetic draw")
  lambda_at_events <- rep(1, nrow(ev$events))
  w <- trigger_pair_weights(ev, pairs, A = 0.3, lambda_at_events, tt, ts)

  # should run without dimension mismatches on an asymmetric x/y grid
  ts2 <- trigger_reestimate_spatial(ev, pairs, w, ec, ts)
  expect_equal(dim(ts2$basevalue), c(length(ts$x), length(ts$y)))
})

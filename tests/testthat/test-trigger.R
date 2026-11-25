make_ready_events2 <- function(n = 50, seed = 11) {
  synth <- make_synthetic_data(n = n, seed = seed)
  ev <- prepare_events(synth$points, synth$region, date_col = "event_date",
                        start_date = "2020-01-01")
  ev <- add_adaptive_bandwidth(ev, k = 3, min_dist = 10, min_bandwidth = 50)
  ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
  ev
}

test_that("trigger_temporal_init returns a kernel integrating to 1", {
  tt <- trigger_temporal_init()
  expect_s3_class(tt, "semihawkes_trigger_temporal")
  expect_equal(simpson(tt$basevalue, 0.05), 1, tolerance = 1e-6)
  expect_equal(tt$fun(-5), 0)
  expect_equal(tt$fun(200), 0)
  expect_true(tt$fun(1) > 0)
})

test_that("trigger_spatial_init returns a kernel integrating to 1 and is symmetric", {
  ts <- trigger_spatial_init()
  expect_s3_class(ts, "semihawkes_trigger_spatial")
  expect_equal(simpson.2D(ts$basevalue, dx = 200, dy = 200), 1, tolerance = 1e-6)
  # isotropic kernel: value at (200,0) should equal value at (0,200) and (-200,0)
  expect_equal(ts$fun(200, 0), ts$fun(0, 200), tolerance = 1e-8)
  expect_equal(ts$fun(200, 0), ts$fun(-200, 0), tolerance = 1e-8)
  # far outside the grid -> 0, not NA
  expect_equal(ts$fun(1e7, 1e7), 0)
})

test_that("trigger_spatial_init handles asymmetric x/y grids correctly (no reshape bug)", {
  ts <- trigger_spatial_init(x_range = 1000, y_range = 2000, x_step = 100, y_step = 250)
  expect_equal(length(ts$x), 21) # -1000 to 1000 by 100
  expect_equal(length(ts$y), 17) # -2000 to 2000 by 250
  expect_equal(dim(ts$basevalue), c(21, 17))
  # value should depend only on (x/scale)^2 + (y/scale)^2 -- check the closed
  # form directly against a couple of grid points, unaffected by the
  # nx != ny mismatch.
  ix <- which(ts$x == 300)
  iy <- which(ts$y == -500)
  raw <- 1 / ((300^2) / (40000^2) + (500^2) / (40000^2) + 0.01)
  raw <- raw / simpson.2D(1 / ((ts$x %o% rep(1, length(ts$y)))^2 / 40000^2 +
                                 (rep(1, length(ts$x)) %o% ts$y)^2 / 40000^2 + 0.01),
                           dx = 100, dy = 250)
  expect_equal(ts$basevalue[ix, iy], raw, tolerance = 1e-8)
})

test_that("the kernel initialisers reject grids Simpson's rule cannot integrate", {
  # seq(-1500, 1500, 200) has 16 points and misses 0: simpson.2D() used to
  # fail with "non-conformable arguments".
  expect_error(trigger_spatial_init(x_range = 1500, x_step = 200),
               "`x_range / x_step` must be a whole number.*give 7.5")
  expect_error(trigger_spatial_init(x_range = 1000, y_range = 1100),
               "`y_range / y_step` must be a whole number")
  # seq(0, 10, 0.4) has 26 points: the kernel used to integrate to 1.016,
  # with only a recycling warning.
  expect_error(trigger_temporal_init(cutoff = 10, step = 0.4),
               "`cutoff / step` must be an even whole number.*give 25")
  expect_error(trigger_temporal_init(cutoff = 10, step = 0.3), "give 33.3333")
  expect_error(trigger_temporal_init(step = 0), "single positive number")
  expect_error(trigger_spatial_init(x_step = -200), "single positive number")

  # Grids that fit still integrate to 1 and, for space, include the origin.
  tt <- trigger_temporal_init(cutoff = 10, step = 0.5)
  expect_equal(simpson(tt$basevalue, 0.5), 1)
  ts <- trigger_spatial_init(x_range = 1400, x_step = 200)
  expect_length(ts$x, 15)
  expect_true(0 %in% ts$x)
  expect_equal(simpson.2D(ts$basevalue, dx = 200, dy = 200), 1)
})

test_that("trigger_evaluate with exclude_self = TRUE zeroes the diagonal contribution", {
  ev <- make_ready_events2(n = 20)
  bg <- background_spatial(ev, cell_size = 150)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()

  res_excl <- trigger_evaluate(ev, bg, tt, ts, exclude_self = TRUE)
  res_incl <- trigger_evaluate(ev, bg, tt, ts, exclude_self = FALSE)

  # excluding self should give strictly smaller (or equal, if the self term
  # happened to be 0) at-events contributions everywhere
  expect_true(all(res_excl$at_events <= res_incl$at_events + 1e-12))
  # the difference at each event i should equal g_T(0) * g_S(0,0)
  self_term <- tt$fun(0) * ts$fun(0, 0)
  expect_equal(res_incl$at_events - res_excl$at_events, rep(self_term, nrow(ev$events)),
               tolerance = 1e-8)
})

test_that("trigger_evaluate matches a brute-force double loop", {
  ev <- make_ready_events2(n = 15)
  bg <- background_spatial(ev, cell_size = 200)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()

  res <- trigger_evaluate(ev, bg, tt, ts, exclude_self = TRUE)

  n <- nrow(ev$events)
  brute <- numeric(n)
  for (j in seq_len(n)) {
    total <- 0
    for (i in seq_len(n)) {
      if (i == j) next
      total <- total + tt$fun(ev$events$days[j] - ev$events$days[i]) *
        ts$fun(ev$events$coorx[j] - ev$events$coorx[i],
               ev$events$coory[j] - ev$events$coory[i])
    }
    brute[j] <- total
  }
  expect_equal(res$at_events, brute, tolerance = 1e-8)
})

test_that("trigger_evaluate's at_all matches integrating over every time mark and cell", {
  ev <- make_ready_events2(n = 30)
  bg <- background_spatial(ev, cell_size = 50)
  tt <- trigger_temporal_init()                          # 90 days, in a one-year window
  ts <- trigger_spatial_init(x_range = 300, x_step = 50) # 300 m, in a 1 km square

  res <- trigger_evaluate(ev, bg, tt, ts)

  d <- ev$events
  grid_x <- bg$x %o% rep(1, length(bg$y))
  grid_y <- rep(1, length(bg$x)) %o% bg$y
  cell_area <- as.numeric(sf::st_area(ev$region)) / sum(bg$marks > 0)
  brute <- 0
  for (i in seq_len(nrow(d))) {
    brute <- brute + mean(tt$fun(ev$time.marks - d$days[i])) * ev$TT *
      sum((bg$marks > 0) * ts$fun(grid_x - d$coorx[i], grid_y - d$coory[i])) * cell_area
  }
  expect_equal(res$at_all, brute, tolerance = 1e-12)
})

test_that("trigger_evaluate gives the same result, to the last bit, with any number of workers", {
  # The events are split into the same chunks whatever the number of workers
  # (four here), and the chunks' sums are added up in order.
  ev <- make_ready_events2(n = 200)
  bg <- background_spatial(ev, cell_size = 200)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()
  expect_length(.fixed_chunks(200), 4)

  res1 <- trigger_evaluate(ev, bg, tt, ts, nprocs = 1)
  res2 <- trigger_evaluate(ev, bg, tt, ts, nprocs = 2)
  expect_identical(res2, res1)

  # A cluster passed in is used and left running, with nothing left behind
  # in its workers.
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  expect_identical(trigger_evaluate(ev, bg, tt, ts, nprocs = cl), res1)
  expect_identical(trigger_evaluate(ev, bg, tt, ts, nprocs = cl), res1)
  expect_equal(unlist(parallel::clusterEvalQ(cl, ls(globalenv(), all.names = TRUE))),
               character(0))
})

test_that("trigger_evaluate validates its inputs' classes", {
  ev <- make_ready_events2(n = 10)
  bg <- background_spatial(ev, cell_size = 200)
  tt <- trigger_temporal_init()
  ts <- trigger_spatial_init()

  expect_error(trigger_evaluate("nope", bg, tt, ts), "semihawkes_events")
  expect_error(trigger_evaluate(ev, "nope", tt, ts), "semihawkes_spatial_base")
  expect_error(trigger_evaluate(ev, bg, "nope", ts), "semihawkes_trigger_temporal")
  expect_error(trigger_evaluate(ev, bg, tt, "nope"), "semihawkes_trigger_spatial")
})

# Sets environment variables or options for the duration of `code` only.
with_envvars <- function(vars, code) {
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(old[nm]))
    }
  })
  do.call(Sys.setenv, as.list(vars))
  force(code)
}
with_opts <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old))
  force(code)
}

test_that("nprocs must be a whole number of workers, \"auto\" or a cluster", {
  for (bad in list("two", 0, 1.5, c(1, 2), NA_real_, -1)) {
    expect_error(.check_nprocs(bad), "whole number of workers")
  }
  expect_silent(.check_nprocs(1))
  expect_silent(.check_nprocs(3L))
  expect_silent(.check_nprocs("auto"))
  ev <- make_ready_events2(n = 10)
  bg <- background_spatial(ev, cell_size = 200)
  expect_error(trigger_evaluate(ev, bg, trigger_temporal_init(), trigger_spatial_init(),
                                nprocs = "two"), "whole number of workers")
})

test_that("a number of workers above the machine's cores is reduced, with a warning", {
  cores <- parallel::detectCores()
  skip_if(is.na(cores))
  expect_warning(n <- .resolve_workers(cores + 5), "more than the")
  expect_equal(n, cores)
  expect_equal(.resolve_workers(1), 1L)
})

test_that("\"auto\" chooses workers to suit the machine", {
  cores <- parallel::detectCores(logical = FALSE)
  skip_if(is.na(cores))
  no_check <- c("_R_CHECK_LIMIT_CORES_" = "", "_R_CHECK_PACKAGE_NAME_" = "")

  with_envvars(no_check, with_opts(list(mc.cores = NULL), {
    # One per physical core but one, though 2 on a 2-core machine.
    expect_equal(.auto_workers(), max(min(2L, cores), cores - 1L))
    # Serial for small data, where starting workers costs more than they save.
    expect_equal(.resolve_workers("auto", n = 999), 1L)
    expect_equal(.resolve_workers("auto", n = 1000), .auto_workers())
  }))
  # At most 2 while R CMD check runs, as CRAN requires.
  with_envvars(c("_R_CHECK_LIMIT_CORES_" = "TRUE", "_R_CHECK_PACKAGE_NAME_" = ""),
               expect_lte(.auto_workers(), 2L))
  with_envvars(c("_R_CHECK_LIMIT_CORES_" = "", "_R_CHECK_PACKAGE_NAME_" = "spHawkes"),
               expect_lte(.auto_workers(), 2L))
  # No more than a user's `mc.cores` option.
  with_envvars(no_check, with_opts(list(mc.cores = 1), expect_equal(.auto_workers(), 1L)))
})

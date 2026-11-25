#' Build the (i, j) triggering-pair index for re-estimation
#'
#' Enumerates ordered event pairs `(i, j)` where `j` could plausibly have
#' triggered `i`: `j` happened strictly before `i`, within `window_days`, and
#' within `window_dist` in both spatial dimensions.
#'
#' @param events A `semihawkes_events` object (from [prepare_events()]).
#' @param window_days Maximum number of days after `j` that `i` can occur;
#'   default `90` (matches the original script's actual filter, `<= 90.0`).
#' @param window_dist Maximum absolute distance (in both x and y) between
#'   `i` and `j`, in the linear units of the working CRS (metres for a UTM
#'   projection); default `1000`, i.e. 1 km, the same as the default reach of
#'   the spatial triggering kernel (see *Choosing the triggering reach* in
#'   [fit_semihawkes()]).
#'
#' @return A data frame with integer columns `i`, `j` (row indices into
#'   `events$events`), one row per qualifying pair.
#'
#' @section Scalability:
#' Once the events are sorted by day, the events that `j` could have
#' triggered (those after it and within `window_days`) form one run, found
#' with [findInterval()]; only those are checked against `window_dist`. Time
#' and memory therefore grow with the number of candidate pairs rather than
#' with `n^2`, as they did when every pair was compared. Pairs are returned
#' in the same order as before: by `j`, then by `i`.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' pairs <- trigger_pairs(ev, window_days = 30, window_dist = 1000)
#' nrow(pairs)
#' head(pairs)
#' @export
trigger_pairs <- function(events, window_days = 90, window_dist = 1000) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  ev <- events$events
  days <- ev$days
  x <- ev$coorx
  y <- ev$coory

  # In time order, the events j could have triggered run from the first one
  # after day days[j] to the last one on or before days[j] + window_days.
  by_day <- order(days)
  sorted_days <- days[by_day]
  first <- findInterval(days, sorted_days) + 1L
  last <- findInterval(days + window_days, sorted_days)

  found <- lapply(seq_along(days), function(j) {
    if (first[j] > last[j]) return(integer(0))
    i <- by_day[first[j]:last[j]]
    i[abs(x[i] - x[j]) <= window_dist & abs(y[i] - y[j]) <= window_dist]
  })
  i <- as.integer(unlist(found, use.names = FALSE))
  j <- rep(seq_along(days), lengths(found))

  ord <- order(j, i) # by j, then i, as before
  data.frame(i = i[ord], j = j[ord])
}

#' Cache per-event spatial edge-correction marks
#'
#' For each event `i`, shifts the triggering kernel's spatial grid
#' (`spatial_trigger$x`/`spatial_trigger$y`) by that event's coordinates and
#' flags which shifted grid cells fall inside the study region, caching the
#' result to disk (one `.rds` file per event). [trigger_edge_correction()]
#' reads these back to compute each event's spatial edge-correction factor
#' and the pooled spatial "repetance" surface.
#'
#' A grid cell counts as inside if it intersects the region, as
#' [sf::st_intersects()] decides, so cells on the boundary count as inside.
#' To keep this fast for many events, the shifted grids of up to 200 events
#' are tested at once, mostly with spatstat's point-in-polygon test; cells
#' within a hair's breadth of the boundary (a thousandth of the grid step)
#' are tested with [sf::st_intersects()] itself, so the marks are the same
#' as testing every cell with sf.
#'
#' @param events A `semihawkes_events` object (from [prepare_events()]).
#' @param spatial_trigger A `semihawkes_trigger_spatial` object (from
#'   [trigger_spatial_init()] or [trigger_reestimate_spatial()]) -- its
#'   `x`/`y` grid defines the marks that get cached.
#' @param cache_dir Directory to write/read `.rds` mark files; created if it
#'   doesn't exist. Required, no default.
#' @param force If `TRUE`, recompute and overwrite marks that are already
#'   cached; default `FALSE` (skip events with an existing file, matching
#'   the original script's `if (!file.exists(fn))` guard).
#' @param nprocs Parallel workers, as for [trigger_evaluate()]: a whole
#'   number (the default, `1`, runs serially), `"auto"`, or a cluster from
#'   [parallel::makeCluster()], which is used and left running.
#'
#' @return `cache_dir`, invisibly.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' ts <- trigger_spatial_init(x_range = 2000, x_step = 250)
#' cache_dir <- tempfile("spHawkes_marks_")
#' trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
#' length(list.files(cache_dir))  # one file per event
#' @export
trigger_spatial_marks <- function(events, spatial_trigger, cache_dir,
                                  force = FALSE, nprocs = 1) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  if (!inherits(spatial_trigger, "semihawkes_trigger_spatial")) {
    stop("`spatial_trigger` must be a `semihawkes_trigger_spatial` object.",
         call. = FALSE)
  }
  if (missing(cache_dir) || is.null(cache_dir)) {
    stop("`cache_dir` is required (no default) -- pass a directory to cache ",
         "the per-event spatial marks in.", call. = FALSE)
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  ev <- events$events
  n <- nrow(ev)
  region <- events$region
  nx <- length(spatial_trigger$x)
  ny <- length(spatial_trigger$y)
  Xv <- as.vector(spatial_trigger$x %o% rep(1, ny))
  Yv <- as.vector(rep(1, nx) %o% spatial_trigger$y)
  ng <- length(Xv)

  # Cells within `eps` of the region's boundary are tested with
  # sf::st_intersects() itself; see .region_point_test().
  eps <- 1e-3 * min(spatial_trigger$x[2] - spatial_trigger$x[1],
                    spatial_trigger$y[2] - spatial_trigger$y[1])
  in_region <- .region_point_test(region, eps)

  files <- file.path(cache_dir, sprintf("mark_%06d.rds", seq_len(n)))

  compute_chunk <- function(idx) {
    todo <- if (force) idx else idx[!file.exists(files[idx])]
    # Up to 200 events' shifted grids per test, to keep memory modest.
    for (block in split(todo, ceiling(seq_along(todo) / 200))) {
      inside <- in_region(rep(Xv, length(block)) + rep(ev$coorx[block], each = ng),
                          rep(Yv, length(block)) + rep(ev$coory[block], each = ng))
      for (k in seq_along(block)) {
        saveRDS(matrix(inside[(k - 1) * ng + seq_len(ng)], nrow = nx, ncol = ny),
                files[block[k]])
      }
    }
    invisible(NULL)
  }

  # `compute_chunk` takes the variables it uses to the workers with it, in its
  # environment, so nothing needs exporting to them; each worker writes its
  # events' files itself.
  workers <- .open_workers(nprocs, n)
  if (is.null(workers$cl)) {
    compute_chunk(seq_len(n))
  } else {
    if (workers$owned) on.exit(parallel::stopCluster(workers$cl), add = TRUE)
    parallel::parLapply(workers$cl, parallel::splitIndices(n, length(workers$cl)),
                        compute_chunk)
  }

  invisible(cache_dir)
}

#' @keywords internal
#' @noRd
.read_spatial_mark <- function(cache_dir, i) {
  fn <- file.path(cache_dir, sprintf("mark_%06d.rds", i))
  if (!file.exists(fn)) {
    stop("Spatial mark for event ", i, " not found in `cache_dir`. Call ",
         "`trigger_spatial_marks()` first.", call. = FALSE)
  }
  readRDS(fn)
}

#' Read every event's spatial mark from `cache_dir` at once
#'
#' @return A list of `n` logical matrices, the `i`-th being event `i`'s mark.
#' @keywords internal
#' @noRd
.read_spatial_marks <- function(cache_dir, n) {
  lapply(seq_len(n), function(i) .read_spatial_mark(cache_dir, i))
}

#' Per-event temporal/spatial edge-correction factors for trigger re-estimation
#'
#' Computes, for every event, how much of its triggering kernels' mass
#' would have been observable given the study period's end date and the
#' study region's boundary -- the correction factors the original script
#' calls `excite.temporal.edge.correction` / `excite.spatial.edge.correction`
#' -- plus the pooled "repetance" surfaces (`temporal.repetance` /
#' `spatial.repetance` in the original) used to correct the re-estimated
#' kernels for how many events *could* have produced an observation at each
#' offset, given the same boundaries.
#'
#' Requires [trigger_spatial_marks()] to have already populated `cache_dir`
#' for every event in `events`.
#'
#' @param events A `semihawkes_events` object.
#' @param temporal_trigger A `semihawkes_trigger_temporal` object.
#' @param spatial_trigger A `semihawkes_trigger_spatial` object.
#' @param cache_dir Directory previously populated by
#'   [trigger_spatial_marks()] for the same `events`/`spatial_trigger`.
#' @param marks Optional: the marks in `cache_dir`, already read into a list
#'   with one logical matrix per event (as [fit_semihawkes()] does once per
#'   fit, so they are not read from disk at every iteration). `NULL` (the
#'   default) reads them from `cache_dir`.
#'
#' @return A list of class `"semihawkes_trigger_edge_correction"` with
#'   `temporal_edge_correction`, `spatial_edge_correction` (numeric vectors,
#'   one value per event), `temporal_repetance` (numeric vector on
#'   `temporal_trigger$base`), `spatial_repetance` (matrix on
#'   `spatial_trigger$x`/`y`), and `temporal_repetance_fun`/
#'   `spatial_repetance_fun` (ready-to-use evaluators).
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' tt <- trigger_temporal_init()
#' ts <- trigger_spatial_init(x_range = 2000, x_step = 250)
#' cache_dir <- tempfile("spHawkes_marks_")
#' trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
#' ec <- trigger_edge_correction(ev, tt, ts, cache_dir)
#' ec
#' @export
trigger_edge_correction <- function(events, temporal_trigger, spatial_trigger,
                                    cache_dir, marks = NULL) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  if (!inherits(temporal_trigger, "semihawkes_trigger_temporal")) {
    stop("`temporal_trigger` must be a `semihawkes_trigger_temporal` object.",
         call. = FALSE)
  }
  if (!inherits(spatial_trigger, "semihawkes_trigger_spatial")) {
    stop("`spatial_trigger` must be a `semihawkes_trigger_spatial` object.",
         call. = FALSE)
  }

  ev <- events$events
  n <- nrow(ev)
  TT <- events$TT
  step <- temporal_trigger$base[2] - temporal_trigger$base[1]
  nx <- length(spatial_trigger$x)
  ny <- length(spatial_trigger$y)
  dx <- spatial_trigger$x[2] - spatial_trigger$x[1]
  dy <- spatial_trigger$y[2] - spatial_trigger$y[1]

  # Small offset avoids landing exactly on a knot of `temporal_trigger$fun`,
  # matching the original script's `+0.6e-5`. The kernel is zero beyond the
  # last point of its grid, so the sum stops there rather than at TT.
  reach <- max(temporal_trigger$base) + step
  temporal_edge_correction <- vapply(seq_len(n), function(i) {
    sum(temporal_trigger$fun(seq(0, min(TT - ev$days[i], reach), step) + 6e-6))
  }, numeric(1))

  temporal_repetance <- vapply(temporal_trigger$base, function(b) {
    sum(b < TT - ev$days)
  }, numeric(1))
  temporal_repetance_fun <- stats::approxfun(temporal_trigger$base, temporal_repetance,
                                             yleft = 1, yright = 1)

  if (is.null(marks)) {
    marks <- .read_spatial_marks(cache_dir, n)
  } else if (length(marks) != n || !identical(dim(marks[[1]]), c(nx, ny))) {
    stop("`marks` does not match `events` and `spatial_trigger`; leave it as ",
         "NULL to read the marks from `cache_dir`.", call. = FALSE)
  }

  spatial_repetance <- matrix(0, nrow = nx, ncol = ny)
  spatial_edge_correction <- numeric(n)
  for (i in seq_len(n)) {
    mark <- marks[[i]]
    spatial_repetance <- spatial_repetance + mark
    spatial_edge_correction[i] <- simpson.2D(mark * spatial_trigger$basevalue, dx = dx, dy = dy)
  }
  spatial_repetance_fun <- .interp_zero(spatial_trigger$x, spatial_trigger$y, spatial_repetance)

  out <- list(
    temporal_edge_correction = temporal_edge_correction,
    spatial_edge_correction = spatial_edge_correction,
    temporal_repetance = temporal_repetance,
    spatial_repetance = spatial_repetance,
    temporal_repetance_fun = temporal_repetance_fun,
    spatial_repetance_fun = spatial_repetance_fun
  )
  class(out) <- "semihawkes_trigger_edge_correction"
  out
}

#' @export
print.semihawkes_trigger_edge_correction <- function(x, ...) {
  cat("<semihawkes_trigger_edge_correction>\n")
  cat(sprintf("  %d events\n", length(x$temporal_edge_correction)))
  invisible(x)
}

#' EM responsibility weights for triggering pairs
#'
#' Computes `A * g_T(t_i - t_j) * g_S(x_i - x_j, y_i - y_j) / lambda(i)` for
#' every pair in `pairs` -- the numerator of the probability that event `j`
#' triggered event `i`, i.e. the original script's `excite.wghs`. This is
#' shared between the temporal and spatial re-estimation steps
#' ([trigger_reestimate_temporal()] / [trigger_reestimate_spatial()]), so it
#' is computed once here rather than duplicated inside each.
#'
#' @param events A `semihawkes_events` object.
#' @param pairs A pair index from [trigger_pairs()].
#' @param A The current branching-ratio parameter estimate (scalar).
#' @param lambda_at_events The current total-intensity estimate at every
#'   event (numeric vector, `length(lambda_at_events) == nrow(events$events)`
#'   -- e.g. `background$bg_at_events + A * trigger_evaluate(...)$at_events`).
#' @param temporal_trigger,spatial_trigger The *current* triggering-kernel
#'   estimates (from the previous EM iteration) used to evaluate the
#'   responsibility weights.
#'
#' @return A numeric vector, one weight per row of `pairs`.
#' @examples
#' \donttest{
#' # A quick fit to one year of the bundled simulated data: a coarse grid, a
#' # coarse triggering-kernel grid and two EM iterations keep it to a few
#' # seconds. Real analyses need finer settings and more iterations.
#' cache_dir <- tempfile("spHawkes_marks_")
#' fit <- fit_semihawkes(
#'   sim_events, sim_region,
#'   date_col  = "date",
#'   end_date  = "2022-01-01",
#'   cell_size = 1000,
#'   cache_dir = cache_dir,
#'   trigger_spatial_args = list(x_range = 2000, x_step = 250),
#'   trigger_window_dist  = 2000,
#'   n_iter  = 2,
#'   verbose = FALSE,
#'   seed    = 1
#' )
#' pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
#' w <- trigger_pair_weights(fit$events, pairs, fit$A, fit$lambda_at_events,
#'                           fit$temporal_trigger, fit$spatial_trigger)
#' summary(w)
#' }
#' @export
trigger_pair_weights <- function(events, pairs, A, lambda_at_events,
                                 temporal_trigger, spatial_trigger) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  ev <- events$events
  A * temporal_trigger$fun(ev$days[pairs$i] - ev$days[pairs$j]) *
    spatial_trigger$fun(ev$coorx[pairs$i] - ev$coorx[pairs$j],
                        ev$coory[pairs$i] - ev$coory[pairs$j]) /
    lambda_at_events[pairs$i]
}

#' Re-estimate the temporal triggering kernel from EM responsibility weights
#'
#' Generalizes the original script's Loop-1 temporal triggering
#' re-estimation: bins triggering-pair offsets into a weighted histogram
#' (edge-corrected by [trigger_edge_correction()]'s output), smooths with
#' [ker.smooth.conv()], and re-normalizes -- returning a
#' `"semihawkes_trigger_temporal"` object, a drop-in replacement for
#' [trigger_temporal_init()]'s output (same `base` grid, refreshed
#' `basevalue`/`fun`).
#'
#' @param events A `semihawkes_events` object.
#' @param pairs A pair index from [trigger_pairs()].
#' @param weights Responsibility weights from [trigger_pair_weights()]
#'   (same length/order as `pairs`).
#' @param edge_correction A `semihawkes_trigger_edge_correction` object from
#'   [trigger_edge_correction()].
#' @param temporal_trigger The current `semihawkes_trigger_temporal` object
#'   (only its `base` grid is reused; `basevalue`/`fun` are replaced).
#' @param bandwidth Gaussian smoothing bandwidth passed to
#'   [ker.smooth.conv()]; default `5` (matches the original script).
#'
#' @return A `"semihawkes_trigger_temporal"` object.
#' @examples
#' \donttest{
#' # A quick fit to one year of the bundled simulated data: a coarse grid, a
#' # coarse triggering-kernel grid and two EM iterations keep it to a few
#' # seconds. Real analyses need finer settings and more iterations.
#' cache_dir <- tempfile("spHawkes_marks_")
#' fit <- fit_semihawkes(
#'   sim_events, sim_region,
#'   date_col  = "date",
#'   end_date  = "2022-01-01",
#'   cell_size = 1000,
#'   cache_dir = cache_dir,
#'   trigger_spatial_args = list(x_range = 2000, x_step = 250),
#'   trigger_window_dist  = 2000,
#'   n_iter  = 2,
#'   verbose = FALSE,
#'   seed    = 1
#' )
#' pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
#' w <- trigger_pair_weights(fit$events, pairs, fit$A, fit$lambda_at_events,
#'                           fit$temporal_trigger, fit$spatial_trigger)
#' ec <- trigger_edge_correction(fit$events, fit$temporal_trigger,
#'                               fit$spatial_trigger, cache_dir)
#' g_new <- trigger_reestimate_temporal(fit$events, pairs, w, ec,
#'                                      fit$temporal_trigger)
#' plot(g_new$base, g_new$basevalue, type = "l", xlim = c(0, 30),
#'      xlab = "Days since event", ylab = "Density")
#' }
#' @export
trigger_reestimate_temporal <- function(events, pairs, weights, edge_correction,
                                        temporal_trigger, bandwidth = 5) {
  if (!inherits(edge_correction, "semihawkes_trigger_edge_correction")) {
    stop("`edge_correction` must come from `trigger_edge_correction()`.", call. = FALSE)
  }
  ev <- events$events
  offsets <- ev$days[pairs$i] - ev$days[pairs$j]
  denom <- edge_correction$temporal_edge_correction[pairs$j] *
    edge_correction$temporal_repetance_fun(offsets)

  series <- hist.weighted(offsets, .ratio_or_zero(weights, denom),
                          breaks = temporal_trigger$base)
  basevalue <- ker.smooth.conv(series$mids, series$density, bandwidth = bandwidth)

  step <- temporal_trigger$base[2] - temporal_trigger$base[1]
  basevalue <- .normalise_kernel(basevalue, simpson(basevalue, step), "temporal")
  fun <- stats::approxfun(temporal_trigger$base + 1e-12, basevalue, yleft = 0, yright = 0)

  out <- list(base = temporal_trigger$base, basevalue = basevalue, fun = fun)
  class(out) <- "semihawkes_trigger_temporal"
  out
}

#' Re-estimate the spatial triggering kernel from EM responsibility weights
#'
#' Spatial triggering
#' re-estimation: bins triggering-pair offsets into a weighted 2D histogram
#' (edge-corrected by [trigger_edge_correction()]'s output), smooths with
#' [ker.smooth.2D.fft()], masks to a maximum radius, and re-normalizes --
#' returning a `"semihawkes_trigger_spatial"` object, a drop-in replacement
#' for [trigger_spatial_init()]'s output.
#'
#'
#' @param events A `semihawkes_events` object.
#' @param pairs A pair index from [trigger_pairs()].
#' @param weights Responsibility weights from [trigger_pair_weights()].
#' @param edge_correction A `semihawkes_trigger_edge_correction` object.
#' @param spatial_trigger The current `semihawkes_trigger_spatial` object
#'   (its `x`/`y` grid is reused; `basevalue`/`fun` are replaced).
#' @param x_bandwidth,y_bandwidth Gaussian smoothing bandwidths passed to
#'   [ker.smooth.2D.fft()], in the linear units of the working CRS (metres for a UTM projection);
#'   default `250` for both, i.e. 250 m (matches the original).
#' @param mark_radius Radius beyond which the re-estimated kernel is forced
#'   to `0`, in the linear units of the working CRS (metres for a UTM
#'   projection). The default, `NULL`, uses 1.05 times the reach of the
#'   kernel's grid (the smaller of its largest `x` and `y` offsets), which
#'   cuts the kernel to a disc just wider than the grid's half-width: 1050 for
#'   the default +/-1 km grid.
#'
#' @return A `"semihawkes_trigger_spatial"` object.
#'
#' @section Offsets no pair could produce:
#' At grid cells where no event could have produced that spatial offset
#' while staying inside the region (`repetance == 0`), the weighted density
#' is zero too, and the ratio of the two is set to `0` rather than `NaN`.
#' Left as `0 / 0`, a single such cell would spoil the whole kernel through
#' the FFT smoother.
#' @examples
#' \donttest{
#' # A quick fit to one year of the bundled simulated data: a coarse grid, a
#' # coarse triggering-kernel grid and two EM iterations keep it to a few
#' # seconds. Real analyses need finer settings and more iterations.
#' cache_dir <- tempfile("spHawkes_marks_")
#' fit <- fit_semihawkes(
#'   sim_events, sim_region,
#'   date_col  = "date",
#'   end_date  = "2022-01-01",
#'   cell_size = 1000,
#'   cache_dir = cache_dir,
#'   trigger_spatial_args = list(x_range = 2000, x_step = 250),
#'   trigger_window_dist  = 2000,
#'   n_iter  = 2,
#'   verbose = FALSE,
#'   seed    = 1
#' )
#' pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
#' w <- trigger_pair_weights(fit$events, pairs, fit$A, fit$lambda_at_events,
#'                           fit$temporal_trigger, fit$spatial_trigger)
#' ec <- trigger_edge_correction(fit$events, fit$temporal_trigger,
#'                               fit$spatial_trigger, cache_dir)
#' h_new <- trigger_reestimate_spatial(fit$events, pairs, w, ec,
#'                                     fit$spatial_trigger)
#' image(h_new$x, h_new$y, h_new$basevalue, asp = 1,
#'       xlab = "x offset (m)", ylab = "y offset (m)")
#' }
#' @export
trigger_reestimate_spatial <- function(events, pairs, weights, edge_correction,
                                       spatial_trigger, x_bandwidth = 250,
                                       y_bandwidth = 250, mark_radius = NULL) {
  if (!inherits(edge_correction, "semihawkes_trigger_edge_correction")) {
    stop("`edge_correction` must come from `trigger_edge_correction()`.", call. = FALSE)
  }
  if (is.null(mark_radius)) {
    mark_radius <- 1.05 * min(max(spatial_trigger$x), max(spatial_trigger$y))
  } else if (!is.numeric(mark_radius) || length(mark_radius) != 1L ||
             !is.finite(mark_radius) || mark_radius <= 0) {
    stop("`mark_radius` must be a single positive number, or `NULL`.", call. = FALSE)
  }
  ev <- events$events
  dx_ij <- ev$coorx[pairs$i] - ev$coorx[pairs$j]
  dy_ij <- ev$coory[pairs$i] - ev$coory[pairs$j]
  denom <- edge_correction$spatial_edge_correction[pairs$j]

  series <- hist.weighted.2D(dx_ij, dy_ij, .ratio_or_zero(weights, denom),
                             x.breaks = spatial_trigger$x, y.breaks = spatial_trigger$y)

  # Fixed meshgrid (see @section above): y.mids for the y coordinate, not
  # x.mids again.
  x_mesh <- series$x.mids %o% rep(1, length(series$y.mids))
  y_mesh <- rep(1, length(series$x.mids)) %o% series$y.mids
  repetance_at_offsets <- edge_correction$spatial_repetance_fun(x_mesh, y_mesh)

  # Offsets no pair in this region could realise have zero repetance; the
  # density there is zero too, so the ratio is zero rather than `NaN`. Left as
  # `0 / 0` a single cell poisons the whole kernel through the FFT smoother.
  basevalue <- ker.smooth.2D.fft(series$x.mids, series$y.mids,
                                 .ratio_or_zero(series$density, repetance_at_offsets),
                                 x.bandwidth = x_bandwidth, y.bandwidth = y_bandwidth)

  grid_x <- spatial_trigger$x %o% rep(1, length(spatial_trigger$y))
  grid_y <- rep(1, length(spatial_trigger$x)) %o% spatial_trigger$y
  mask <- (grid_x^2 + grid_y^2) < mark_radius^2
  basevalue <- basevalue * mask

  dx <- spatial_trigger$x[2] - spatial_trigger$x[1]
  dy <- spatial_trigger$y[2] - spatial_trigger$y[1]
  basevalue <- .normalise_kernel(basevalue, simpson.2D(basevalue, dx = dx, dy = dy),
                                 "spatial")
  fun <- .interp_zero(spatial_trigger$x, spatial_trigger$y, basevalue)

  out <- list(x = spatial_trigger$x, y = spatial_trigger$y, basevalue = basevalue, fun = fun)
  class(out) <- "semihawkes_trigger_spatial"
  out
}

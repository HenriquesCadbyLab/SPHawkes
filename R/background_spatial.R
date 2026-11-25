#' Non-parametric spatial background surface
#'
#' Estimates the spatial component of the background intensity on a regular
#' grid covering the bounding box of `events$region`, by summing each
#' event's isotropic Gaussian kernel (using its adaptive `bandwidth` and
#' edge-correction `bg.integral`, see [add_adaptive_bandwidth()] and
#' [add_bg_integral()]) evaluated at every grid cell, then mean-normalizing
#' over the cells that fall inside `region`.
#'
#' As with [temporal_base_weekly()], [temporal_base_monthly()] and
#' [temporal_base_trend()], `resp_weight = NULL` reproduces the EM
#' initial-guess step (every event contributes with weight 1), and supplying
#' `resp_weight` (the E-step background-responsibility vector \eqn{\varphi_i},
#' one value per event) reproduces the weighted re-estimation step.
#'
#' @param events A `semihawkes_events` object (from [prepare_events()]),
#'   already passed through [add_adaptive_bandwidth()] and
#'   [add_bg_integral()] (i.e. `events$events$bandwidth` and
#'   `events$events$bg.integral` must exist).
#' @param cell_size Grid cell size, in the same linear units as the working
#'   CRS (e.g. metres for a UTM projection). There is deliberately no
#'   default: pick a value appropriate to your CRS units and the spatial
#'   scale of your data (the original script used `150` for a dataset
#'   projected in metres).
#' @param resp_weight Optional numeric vector, one value per row of
#'   `events$events`, giving the EM E-step background-responsibility weight
#'   for each event. If `NULL` (the default), every event gets weight `1`
#'   (the initial-guess behaviour).
#' @param chunk_size Number of events processed per vectorized chunk, to cap
#'   memory use for large datasets and/or fine grids (each chunk builds a
#'   `n_grid_cells x chunk_size` matrix internally). Default `200`.
#'
#' @return A list of class `"semihawkes_spatial_base"` with elements:
#'   \describe{
#'     \item{x, y}{The grid's coordinate sequences.}
#'     \item{basevalue}{`length(x) x length(y)` matrix, mean-normalized over
#'       cells with `marks > 0`.}
#'     \item{marks}{`length(x) x length(y)` integer matrix, `1` where the
#'       cell centre falls inside `events$region`, `0` otherwise.}
#'     \item{cell_size}{The `cell_size` used.}
#'     \item{fun}{A function `fun(x, y)` that bilinearly interpolates
#'       `basevalue` at arbitrary coordinates via
#'       [fields::interp.surface()], returning `NA` outside the grid.}
#'   }
#'
#' @references Zhuang, J. and Mateu, J. (2019). A semiparametric
#'   spatiotemporal Hawkes-type point process model with periodic background
#'   for crime data. *Journal of the Royal Statistical Society Series A*,
#'   182(3), 919-942.
#'
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' bg <- background_spatial(ev, cell_size = 1000)
#' bg
#' @export
background_spatial <- function(events, cell_size, resp_weight = NULL,
                               chunk_size = 200) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  if (missing(cell_size) || is.null(cell_size) || !is.numeric(cell_size) ||
      length(cell_size) != 1 || cell_size <= 0) {
    stop("`cell_size` is required: a single positive number, in the same ",
         "linear units as the working CRS (e.g. metres for a UTM ",
         "projection). There is no default because it depends on your ",
         "data's spatial scale and CRS units.", call. = FALSE)
  }
  if (is.null(events$events$bandwidth)) {
    stop("`events$events$bandwidth` is missing; call `add_adaptive_bandwidth()` ",
         "first.", call. = FALSE)
  }
  if (is.null(events$events$bg.integral)) {
    stop("`events$events$bg.integral` is missing; call `add_bg_integral()` ",
         "first.", call. = FALSE)
  }

  ev <- events$events
  n <- nrow(ev)
  if (!is.null(resp_weight) && length(resp_weight) != n) {
    stop("`resp_weight` must have length `nrow(events$events)` (", n, ").",
         call. = FALSE)
  }

  ## --- 1. Build the grid from region's bounding box ---------------------
  bbox <- sf::st_bbox(events$region)
  grid_x <- .cover_seq(bbox[["xmin"]], bbox[["xmax"]], cell_size)
  grid_y <- .cover_seq(bbox[["ymin"]], bbox[["ymax"]], cell_size)
  nx <- length(grid_x)
  ny <- length(grid_y)

  basex <- grid_x %o% rep(1, ny) # nx x ny
  basey <- rep(1, nx) %o% grid_y # nx x ny

  ## --- 2. Which cells fall inside the region? ----------------------------
  grid_points <- expand.grid(x = grid_x, y = grid_y) # x varies fastest
  grid_sf <- sf::st_as_sf(grid_points, coords = c("x", "y"), crs = events$crs)
  inside <- lengths(sf::st_intersects(grid_sf, events$region)) > 0
  marks <- matrix(as.integer(inside), nrow = nx, ncol = ny)

  ## --- 3. Accumulate each event's edge-corrected Gaussian kernel on the
  ## grid, chunked and vectorized over events (grid dimension is already
  ## vectorized; chunking only bounds memory when many events are summed).
  gx <- as.vector(basex)
  gy <- as.vector(basey)
  ng <- length(gx)

  basevalue <- numeric(ng)
  weight <- if (is.null(resp_weight)) rep(1, n) else resp_weight

  chunk_starts <- seq(1, n, by = chunk_size)
  for (start in chunk_starts) {
    ix <- start:min(start + chunk_size - 1, n)
    nk <- length(ix)

    dx <- outer(gx, ev$coorx[ix], "-") # ng x nk
    dy <- outer(gy, ev$coory[ix], "-") # ng x nk
    bw_rep <- rep(ev$bandwidth[ix], each = ng)

    dens <- stats::dnorm(dx, sd = bw_rep) * stats::dnorm(dy, sd = bw_rep)
    contrib <- dens / rep(ev$bg.integral[ix], each = ng) * rep(weight[ix], each = ng)

    basevalue <- basevalue + rowSums(matrix(contrib, nrow = ng))
  }

  basevalue <- matrix(basevalue, nrow = nx, ncol = ny)
  basevalue <- basevalue / mean(basevalue[marks > 0])

  fun <- (function(grid_x, grid_y, basevalue) {
    function(x, y) {
      fields::interp.surface(
        obj = list(x = grid_x, y = grid_y, z = basevalue),
        loc = cbind(x = c(x), y = c(y))
      )
    }
  })(grid_x, grid_y, basevalue)

  out <- list(x = grid_x, y = grid_y, basevalue = basevalue, marks = marks,
              cell_size = cell_size, fun = fun)
  class(out) <- "semihawkes_spatial_base"
  out
}

#' @export
print.semihawkes_spatial_base <- function(x, ...) {
  cat("<semihawkes_spatial_base>\n")
  cat(sprintf("  grid: %d x %d cells (cell_size = %s), %d inside `region`\n",
              length(x$x), length(x$y), format(x$cell_size), sum(x$marks > 0)))
  cat(sprintf("  basevalue: mean-normalized over in-region cells (mean = %.4f)\n",
              mean(x$basevalue[x$marks > 0])))
  invisible(x)
}

#' Combine temporal and spatial background components into the background
#' intensity shape (without the scalar `mu` parameter)
#'
#' Multiplies the trend, weekly and monthly temporal base values by the
#' spatial background surface, evaluated at the observed events and summarized
#' over the study grid.
#'
#' @param events A `semihawkes_events` object (from [prepare_events()]).
#' @param trend,weekly,monthly `"semihawkes_temporal_base"` objects from
#'   [temporal_base_trend()], [temporal_base_weekly()] and
#'   [temporal_base_monthly()] respectively.
#' @param spatial A `"semihawkes_spatial_base"` object from
#'   [background_spatial()].
#'
#' @return A list with elements:
#'   \describe{
#'     \item{mub_events}{Spatial background value at each event's location.}
#'     \item{bg_at_events}{`trend * weekly * monthly * mub_events` at each
#'       event -- the background intensity shape at the observed events,
#'       i.e. \eqn{\mu(s_i,t_i)/\mu_0}.}
#'     \item{mub_grid}{Spatial background value on the full spatial grid
#'       (same shape as `spatial$basevalue`).}
#'     \item{bg_at_all}{Single number: the space-time integral of the
#'       background intensity shape over the study region and period,
#'       \eqn{\iint_S \int_0^T \mu(s,t)/\mu_0 \, dt \, ds}, used to normalize
#'       the log-likelihood.}
#'     \item{marks_events}{Calendar marks (weekday and month) for each
#'       event, as returned by `.calendar_marks()`; consumed downstream by
#'       `em_step()`.}
#'   }
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' bgint <- background_intensity(
#'   ev,
#'   trend   = temporal_base_trend(ev),
#'   weekly  = temporal_base_weekly(ev),
#'   monthly = temporal_base_monthly(ev),
#'   spatial = background_spatial(ev, cell_size = 1000)
#' )
#' summary(bgint$bg_at_events)
#' bgint$bg_at_all
#' @export
background_intensity <- function(events, trend, weekly, monthly, spatial) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  for (nm in c("trend", "weekly", "monthly")) {
    obj <- get(nm)
    if (!inherits(obj, "semihawkes_temporal_base")) {
      stop("`", nm, "` must be a `semihawkes_temporal_base` object.", call. = FALSE)
    }
  }
  if (!inherits(spatial, "semihawkes_spatial_base")) {
    stop("`spatial` must be a `semihawkes_spatial_base` object (see ",
         "`background_spatial()`).", call. = FALSE)
  }

  ev <- events$events
  cell_area <- .grid_cell_area(spatial, events$region)

  marks_events <- .calendar_marks(ev$days, events$start_date)
  marks_grid   <- .calendar_marks(events$time.marks, events$start_date)

  mub_events <- spatial$fun(ev$coorx, ev$coory)
  if (anyNA(mub_events)) {
    stop(sum(is.na(mub_events)), " of ", length(mub_events), " events fall ",
         "outside the spatial background grid, so their background intensity ",
         "is undefined. `spatial` was probably built from different events or ",
         "a different region; rebuild it with `background_spatial()` on these ",
         "events.", call. = FALSE)
  }
  bg_at_events <- trend$fun(ev$days) * weekly$fun(marks_events$wday) *
    monthly$fun(marks_events$month) * mub_events

  mub_grid <- spatial$basevalue
  bg_at_all <- as.numeric(
    mean(trend$fun(events$time.marks) * weekly$fun(marks_grid$wday) *
           monthly$fun(marks_grid$month)) * events$TT *
      sum(mub_grid * spatial$marks) * cell_area
  )

  list(mub_events = mub_events, bg_at_events = bg_at_events,
       mub_grid = mub_grid, bg_at_all = bg_at_all,
       marks_events = marks_events)
}

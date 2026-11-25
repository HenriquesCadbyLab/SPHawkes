#' Initial (parametric) temporal triggering-function guess
#'
#' Builds the Loop-0 power-law temporal triggering kernel
#' \eqn{g(t) = (\text{offset} + t/\text{scale})^{\text{exponent}}}, normalized
#' to integrate to `1` over `[0, cutoff]` via [simpson()].
#'
#' @param cutoff Maximum time lag considered, in days; default `90`.
#' @param step Grid step for `base`, in days; default `0.05`. `cutoff / step`
#'   must be an even whole number (1800 for the defaults), so that the grid
#'   ends at `cutoff` and has the odd number of points [simpson()] needs.
#' @param offset,scale,exponent Power-law shape parameters:
#'   \eqn{g(t) = (\text{offset} + t/\text{scale})^{\text{exponent}}}. Defaults
#'   (`0.1`, `10`, `-1.03`) match the original script.
#'
#' @return A list of class `"semihawkes_trigger_temporal"` with `base` (the
#'   time-lag grid), `basevalue` (normalized kernel values), and `fun` (a
#'   vectorized evaluator, `0` outside `[0, cutoff]`).
#' @examples
#' g <- trigger_temporal_init()
#' g
#' plot(g$base, g$basevalue, type = "l", log = "y",
#'      xlab = "Days since event", ylab = "Density")
#' @export
trigger_temporal_init <- function(cutoff = 90, step = 0.05, offset = 0.1,
                                  scale = 10, exponent = -1.03) {
  .grid_steps(cutoff, step, "cutoff", "step", even = TRUE,
              grid = "`seq(0, cutoff, step)` ends at `cutoff`")
  base <- seq(0, cutoff, step)
  basevalue <- (offset + base / scale)^exponent
  basevalue <- basevalue / simpson(basevalue, step)
  fun <- stats::approxfun(base + 1e-12, basevalue, yleft = 0, yright = 0)

  out <- list(base = base, basevalue = basevalue, fun = fun)
  class(out) <- "semihawkes_trigger_temporal"
  out
}

#' @export
print.semihawkes_trigger_temporal <- function(x, ...) {
  cat("<semihawkes_trigger_temporal>\n")
  cat(sprintf("  %d points over [%.2f, %.2f] days\n",
              length(x$base), min(x$base), max(x$base)))
  invisible(x)
}

#' Initial (parametric) spatial triggering-function guess
#'
#' Builds the Loop-0 isotropic spatial triggering kernel
#' \eqn{g(x,y) = 1 / (x^2/\text{scale}^2 + y^2/\text{scale}^2 + \epsilon)} on a
#' regular grid, normalized to integrate to `1` via [simpson.2D()].
#'
#' @param x_range,y_range Half-width of the grid along `x`/`y` (grid spans
#'   `[-range, range]`), in the linear units of the working CRS (metres for a
#'   UTM projection); default `1000` for both, i.e. a +/-1 km window. This is
#'   how far from an event the triggering kernel can reach, so choose it with
#'   care: see *Choosing the triggering reach* in [fit_semihawkes()].
#' @param x_step,y_step Grid step along `x`/`y`, in the linear units of the working CRS (metres for a UTM projection); default `200` for both, i.e. 200 m.
#'   `x_range / x_step` and `y_range / y_step` must be whole numbers (5 for
#'   the defaults), so that each grid is centred on 0 and has the odd number
#'   of points [simpson.2D()] needs.
#' @param scale,epsilon Kernel shape parameters:
#'   \eqn{g(x,y) = 1 / (x^2/\text{scale}^2 + y^2/\text{scale}^2 + \epsilon)}.
#'   Defaults (`40000`, `0.01`) match the original script.
#'
#' @return A list of class `"semihawkes_trigger_spatial"` with `x`, `y` (grid
#'   coordinate sequences), `basevalue` (`length(x) x length(y)` matrix,
#'   normalized), and `fun` (a function `fun(x, y)` that bilinearly
#'   interpolates via [fields::interp.surface()], returning `0` -- not `NA`
#'   -- outside the grid, matching the original `excite.spatial.fun`).
#' @examples
#' h <- trigger_spatial_init(x_range = 2000, x_step = 100)
#' h
#' image(h$x, h$y, h$basevalue, asp = 1,
#'       xlab = "x offset (m)", ylab = "y offset (m)")
#' @export
trigger_spatial_init <- function(x_range = 1000, y_range = x_range,
                                 x_step = 200, y_step = x_step,
                                 scale = 40000, epsilon = 0.01) {
  .grid_steps(x_range, x_step, "x_range", "x_step",
              grid = "`seq(-x_range, x_range, x_step)` is centred on 0")
  .grid_steps(y_range, y_step, "y_range", "y_step",
              grid = "`seq(-y_range, y_range, y_step)` is centred on 0")
  base_x <- seq(-x_range, x_range, x_step)
  base_y <- seq(-y_range, y_range, y_step)
  nx <- length(base_x)
  ny <- length(base_y)

  X <- base_x %o% rep(1, ny) # nx x ny
  Y <- rep(1, nx) %o% base_y # nx x ny

  basevalue <- 1 / ((abs(X)^2) / (scale^2) + (abs(Y)^2) / (scale^2) + epsilon)
  basevalue <- basevalue / simpson.2D(basevalue, dx = x_step, dy = y_step)

  fun <- .interp_zero(base_x, base_y, basevalue)

  out <- list(x = base_x, y = base_y, basevalue = basevalue, fun = fun)
  class(out) <- "semihawkes_trigger_spatial"
  out
}

#' @export
print.semihawkes_trigger_spatial <- function(x, ...) {
  cat("<semihawkes_trigger_spatial>\n")
  cat(sprintf("  %d x %d grid, x in [%.0f, %.0f], y in [%.0f, %.0f]\n",
              length(x$x), length(x$y), min(x$x), max(x$x), min(x$y), max(x$y)))
  invisible(x)
}

#' Bilinear grid interpolator that returns `0` (not `NA`) outside the grid
#'
#' @param grid_x,grid_y Grid coordinate sequences.
#' @param z `length(grid_x) x length(grid_y)` value matrix.
#' @return A function `fun(x, y)`.
#' @keywords internal
#' @noRd
.interp_zero <- function(grid_x, grid_y, z) {
  function(x, y) {
    temp <- fields::interp.surface(
      obj = list(x = grid_x, y = grid_y, z = z),
      loc = cbind(x = c(x), y = c(y))
    )
    temp[is.na(temp)] <- 0
    temp
  }
}

#' Evaluate self-exciting (triggering) contributions at events and over the
#' study domain
#'
#' For every event \eqn{i}, adds its triggering contribution
#' \eqn{g_T(t_j - t_i) \, g_S(x_j - x_i, y_j - y_i)} to every other event
#' \eqn{j} (this is `0` whenever \eqn{t_j \le t_i}, because `temporal_trigger`
#' and `spatial_trigger` vanish for non-positive/out-of-range offsets), and
#' separately integrates its contribution over the full spatial grid and
#' `[0, TT]`. The kernels are zero beyond the ends of their own grids, so
#' that integral only evaluates them at the time marks and background-grid
#' cells within their reach of the event; the rest would add exactly zero.
#'
#' @param events A `semihawkes_events` object (from [prepare_events()]).
#' @param background A `semihawkes_spatial_base` object (from
#'   [background_spatial()]) -- only its grid and `marks` are used here, to
#'   define the spatial domain integrated over.
#' @param temporal_trigger A `semihawkes_trigger_temporal` object (from
#'   [trigger_temporal_init()], or a later Loop-`k` re-estimate).
#' @param spatial_trigger A `semihawkes_trigger_spatial` object (from
#'   [trigger_spatial_init()], or a later Loop-`k` re-estimate).
#' @param exclude_self Whether to zero each event's contribution to itself;
#'   default `TRUE`.
#' @param nprocs Parallel workers: a whole number (the default, `1`, runs
#'   serially), `"auto"` to choose from the machine, or a cluster from
#'   [parallel::makeCluster()], which is used and left running. `"auto"`
#'   uses one worker per physical core but one (two on a two-core machine),
#'   at most 2 while `R CMD check` runs, and runs serially for fewer than
#'   1,000 events, where starting the workers costs more than they save. A
#'   number larger than the machine's cores is reduced to the number of
#'   cores, with a warning. The workers are socket clusters, which work on
#'   Windows, macOS and Linux, and the result is the same, to the last bit,
#'   whatever their number. Starting them takes up to a second, so for
#'   repeated calls pass a cluster; [fit_semihawkes()] does this itself.
#'
#' @return A list with `at_events` (numeric vector, one value per event: the
#'   total triggering contribution received from all other events) and
#'   `at_all` (single number: the same contribution integrated over the
#'   whole study region and `[0, TT]`, summed across all events -- i.e.
#'   `mytriggers.at.events.no.A` / `mytriggers.at.all.no.A` in the original
#'   script, without the `A` scaling parameter).
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
#' ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
#' bg <- background_spatial(ev, cell_size = 1000)
#' trig <- trigger_evaluate(ev, bg, trigger_temporal_init(),
#'                          trigger_spatial_init(x_range = 2000, x_step = 250))
#' summary(trig$at_events)
#' trig$at_all
#' @export
trigger_evaluate <- function(events, background, temporal_trigger, spatial_trigger,
                             exclude_self = TRUE, nprocs = 1) {
  if (!inherits(events, "semihawkes_events")) {
    stop("`events` must be a `semihawkes_events` object (see `prepare_events()`).",
         call. = FALSE)
  }
  if (!inherits(background, "semihawkes_spatial_base")) {
    stop("`background` must be a `semihawkes_spatial_base` object (see ",
         "`background_spatial()`).", call. = FALSE)
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

  cell_area <- .grid_cell_area(background, events$region)

  # Each event's triggering integrated over [0, TT]: the mean of the
  # temporal kernel over the time marks, times TT. The kernel is zero before
  # the event and beyond the last point of its grid, so only the marks in
  # between are evaluated. Events on the same day share the value.
  time_marks <- events$time.marks
  reach_t <- max(temporal_trigger$base)
  event_days <- sort(unique(ev$days))
  # First and last time mark within reach of each day, looked up for all
  # days at once: findInterval() checks that `time_marks` is sorted at
  # every call, which took most of the time when called day by day.
  first_mark <- pmax(1L, findInterval(event_days, time_marks))
  last_mark <- findInterval(event_days + reach_t + 1e-6, time_marks)
  temporal_mass <- vapply(seq_along(event_days), function(j) {
    k <- seq.int(first_mark[j], last_mark[j])
    sum(temporal_trigger$fun(time_marks[k] - event_days[j]))
  }, numeric(1)) / length(time_marks) * events$TT
  temporal_part <- temporal_mass[match(ev$days, event_days)]

  # The spatial kernel is zero outside its grid, so only the background
  # cells whose offset from the event falls within that grid are evaluated.
  work <- list(
    n = n, days = ev$days, x = ev$coorx, y = ev$coory,
    temporal_part = temporal_part, exclude_self = isTRUE(exclude_self),
    temporal_fun = temporal_trigger$fun, spatial_fun = spatial_trigger$fun,
    reach_x = range(spatial_trigger$x), reach_y = range(spatial_trigger$y),
    bg_x = background$x, bg_y = background$y,
    in_region = background$marks > 0, cell_area = cell_area
  )

  # The events are split into the same chunks whatever the number of workers,
  # and the chunks' sums are added up in order, so the result does not depend
  # on `nprocs`.
  workers <- .open_workers(nprocs, n)
  if (workers$owned) on.exit(parallel::stopCluster(workers$cl), add = TRUE)
  results <- .run_chunks(workers$cl, .fixed_chunks(n), .trigger_evaluate_chunk, work)

  list(
    at_events = Reduce(`+`, lapply(results, function(r) r$event)),
    at_all = sum(vapply(results, function(r) r$all, numeric(1)))
  )
}

# The triggering of the events `idx` at every event (`event`, a vector over
# all events) and integrated over the study period and region (`all`, a
# number), for trigger_evaluate(). It runs in parallel workers, so it uses
# only its arguments: `work` holds the events, kernels and grid.
.trigger_evaluate_chunk <- function(idx, work) {
  ev_sum <- numeric(work$n)
  all_sum <- 0
  for (i in idx) {
    contrib <- work$temporal_fun(work$days - work$days[i]) *
      work$spatial_fun(work$x - work$x[i], work$y - work$y[i])
    if (work$exclude_self) contrib[i] <- 0
    ev_sum <- ev_sum + contrib

    dx <- work$bg_x - work$x[i]
    dy <- work$bg_y - work$y[i]
    ix <- which(dx >= work$reach_x[1] & dx <= work$reach_x[2])
    iy <- which(dy >= work$reach_y[1] & dy <= work$reach_y[2])
    if (length(ix) > 0 && length(iy) > 0) {
      grid_spatial <- work$spatial_fun(dx[ix] %o% rep(1, length(iy)),
                                       rep(1, length(ix)) %o% dy[iy])
      all_sum <- all_sum + work$temporal_part[i] *
        sum(work$in_region[ix, iy] * grid_spatial) * work$cell_area
    }
  }
  list(event = ev_sum, all = all_sum)
}

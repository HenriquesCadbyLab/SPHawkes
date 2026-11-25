#' Evaluate background and triggering intensity on a spatial grid over time
#'
#' Computes, for each day in the study period (or a user-supplied set of
#' days), the background intensity \eqn{\mu(s,t)} and the triggering
#' intensity \eqn{A \cdot g(t) \cdot f(s)} on a regular spatial grid covering
#' `city`. This is the building block behind exceedance-set and persistence
#' analyses: it turns a fitted `semihawkes_fit` object into daily snapshots
#' of both intensity components that can later be aggregated over periods
#' (yearly, semestral, quarterly, ...) via [integrate_intensity_period()].
#'
#' Unlike `background_intensity()` and `em_step()`, which evaluate the model
#' only at observed event locations (for likelihood/EM purposes), this
#' function evaluates it on an arbitrary regular grid, independent of the
#' grid resolution `background_spatial()` used internally during fitting.
#' You choose the grid resolution via `cell_size`.
#'
#' `fit$background$fun`, `fit$trend$fun`, `fit$weekly$fun`, `fit$monthly$fun`,
#' `fit$temporal_trigger$fun` and `fit$spatial_trigger$fun` are assumed to be
#' the final, converged fit -- i.e. `fit` is the result of a completed
#' `fit_semihawkes()` run (or the last `em_step()` in a manual loop), not an
#' initial Loop-0 guess.
#'
#' @param fit A `semihawkes_fit` object (or equivalently-shaped list) with
#'   `$events`, `$mu`, `$A`, `$trend`, `$weekly`, `$monthly`, `$background`,
#'   `$temporal_trigger`, `$spatial_trigger`.
#' @param city An `sf` polygon giving the study region. Used both to build
#'   the evaluation grid (from its bounding box) and to mask grid cells
#'   outside the region (via `spatstat.geom::inside.owin()`).
#' @param cell_size Grid spacing, in the same spatial units as `city` (e.g.
#'   metres for UTM data). Default `150`.
#' @param days Integer vector of days-since-`start_date` at which to
#'   evaluate the intensity. Defaults to `0:TT`, where
#'   `TT <- max(fit$events$events$days) + 1`.
#' @param t.lag Only triggering events within this many days before `t`
#'   contribute (events further in the past are dropped for speed, on the
#'   assumption their contribution has decayed to ~0). Default `90`.
#' @param s.lag Only triggering events within this spatial distance of the
#'   grid's bounding box contribute. The default, `NULL`, uses the reach of
#'   the fitted spatial triggering kernel, beyond which its contribution is
#'   zero.
#' @param nprocs Number of parallel workers, via `foreach`/`doParallel`.
#'   `nprocs = 1` (default) runs sequentially with no cluster overhead. With
#'   more, the workers are stopped on exit and the sequential `foreach`
#'   backend is registered in their place.
#'
#' @return A list with:
#'   \describe{
#'     \item{`days`}{the days at which intensity was evaluated (matches
#'       the `days` argument).}
#'     \item{`x`, `y`}{grid coordinate vectors.}
#'     \item{`inside`}{logical vector, `TRUE` for grid cells inside `city`.}
#'     \item{`snapshots`}{a list (one element per day) of lists with `$bg`
#'       and `$tr`, each a `length(x)` by `length(y)` matrix of intensity
#'       values. Cells outside `city` are `0`.}
#'   }
#'
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
#' # Evaluate on a handful of days rather than the whole study window
#' grid <- background_trigger_grid(fit, sim_region, cell_size = 500,
#'                                 days = c(0, 100, 200))
#' length(grid$snapshots)
#' dim(grid$snapshots[[1]]$bg)
#' }
#' @export
background_trigger_grid <- function(fit, city, cell_size = 150,
                                    days = NULL, t.lag = 90, s.lag = NULL,
                                    nprocs = 1) {
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
    stop("Package 'spatstat.geom' is required for background_trigger_grid().",
         call. = FALSE)
  }
  if (is.null(s.lag)) {
    s.lag <- max(abs(c(fit$spatial_trigger$x, fit$spatial_trigger$y)))
  }

  bbox <- sf::st_bbox(city)
  x_seq <- .cover_seq(bbox[["xmin"]], bbox[["xmax"]], cell_size)
  y_seq <- .cover_seq(bbox[["ymin"]], bbox[["ymax"]], cell_size)
  grid <- expand.grid(x = x_seq, y = y_seq)

  w <- spatstat.geom::as.owin(city)
  inside <- spatstat.geom::inside.owin(grid$x, grid$y, w = w)
  idx <- which(inside)

  bg_spatial_vals <- fit$background$fun(grid$x, grid$y)

  history <- as.matrix(fit$events$events[, c("coorx", "coory", "days")])

  if (is.null(days)) {
    TT <- max(fit$events$events$days) + 1
    days <- 0:TT
  }

  # Temporal part of the background for every requested day. `weekly$fun` and
  # `monthly$fun` are periodic in calendar time, so they take the day of the
  # week and the calendar month, not days since `start_date` -- passing the
  # day offset directly shifts the weekly pattern unless the study starts on
  # a Sunday, and makes the "monthly" term repeat every 12 days. Computed
  # here, once, so the per-day worker below needs no package internals.
  marks <- .calendar_marks(days, fit$events$start_date)
  bg_temporal <- fit$mu *
    fit$trend$fun(days) *
    fit$weekly$fun(marks$wday) *
    fit$monthly$fun(marks$month)

  n_grid <- nrow(grid)
  nx <- length(x_seq)
  ny <- length(y_seq)

  evaluate_one_day <- function(t) {
    lb_bg <- numeric(n_grid)
    lb_tr <- numeric(n_grid)

    if (length(idx) == 0L) {
      z <- matrix(0, nx, ny)
      return(list(bg = z, tr = z))
    }

    lb_bg[idx] <- bg_temporal[match(t, days)] * bg_spatial_vals[idx]

    past <- history[history[, 3] < t & (t - history[, 3]) <= t.lag, , drop = FALSE]

    if (nrow(past) > 0L) {
      keep <- past[, 1] >= min(grid$x) - s.lag &
        past[, 1] <= max(grid$x) + s.lag &
        past[, 2] >= min(grid$y) - s.lag &
        past[, 2] <= max(grid$y) + s.lag
      past <- past[keep, , drop = FALSE]
    }

    if (nrow(past) > 0L) {
      # Vectorised over both grid cells (idx) and past events (rows of
      # `past`) via outer differences, instead of accumulating one event
      # at a time in a for-loop.
      dt <- t - past[, 3]                                   # length n_past
      temporal_part <- fit$temporal_trigger$fun(dt)          # length n_past

      dx <- outer(grid$x[idx], past[, 1], "-")                # n_idx x n_past
      dy <- outer(grid$y[idx], past[, 2], "-")                # n_idx x n_past

      # spatial_trigger$fun() flattens its output (loses dim()) even when
      # given matrix input -- restore the n_idx x n_past shape explicitly
      # before reducing over events, rather than assuming it comes back
      # as a matrix.
      spatial_part <- fit$spatial_trigger$fun(dx, dy)
      dim(spatial_part) <- dim(dx)

      contrib <- spatial_part %*% temporal_part               # n_idx x 1
      lb_tr[idx] <- fit$A * as.vector(contrib)
    }

    list(
      bg = matrix(lb_bg, nrow = nx, ncol = ny),
      tr = matrix(lb_tr, nrow = nx, ncol = ny)
    )
  }

  if (nprocs > 1L) {
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for nprocs > 1.",
           call. = FALSE)
    }
    cl <- parallel::makeCluster(nprocs)
    doParallel::registerDoParallel(cl)
    # Unregister the workers as well as stopping them. Left registered, the
    # stopped cluster stays the foreach backend, and the next `%dopar%` in
    # the session, from any package, fails with "invalid connection".
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)

    t <- NULL # avoid a spurious R CMD check NOTE for the foreach binding
    snapshots <- foreach::foreach(
      t = days,
      .packages = c("spatstat.geom")
    ) %dopar% {
      evaluate_one_day(t)
    }
  } else {
    snapshots <- lapply(days, evaluate_one_day)
  }

  list(days = days, x = x_seq, y = y_seq, inside = inside, snapshots = snapshots)
}

#' @importFrom foreach %dopar%
NULL

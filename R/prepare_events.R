#' Prepare a spatio-temporal event dataset for semi-parametric Hawkes fitting
#'
#' Validates, cleans and reshapes a set of georeferenced events (e.g. crime
#' records) and a study-region polygon into the flat structure used by the
#' rest of the `spHawkes` EM pipeline: projected `x`/`y` coordinates, a
#' numeric time index in days since `start_date`, and events restricted to
#' `region` and to `[start_date, end_date)`.
#'
#'
#' @param points An `sf` object with `POINT` geometry (one row per event).
#' @param region An `sf`/`sfc` object with `POLYGON`/`MULTIPOLYGON` geometry
#'   describing the study area (e.g. a city or district boundary).
#' @param date_col Name of the column in `points` holding the event date
#'   (character scalar). The column may be `Date`, `POSIXct`/`POSIXlt`, or a
#'   string parseable by [as.Date()]. Date-times are converted to calendar
#'   dates in their own time zone (their `tzone` attribute, or the session's
#'   time zone if they have none), so an event at 20:00 on a Friday in
#'   Chicago stays on that Friday.
#' @param start_date Start of the study period. If `NULL` (default), the
#'   earliest event date is used. Coerced with [as.Date()], or, for a
#'   date-time, converted in its own time zone as for `date_col`.
#' @param end_date Last date for selecting events, exclusive (events with
#'   `date >= end_date` are dropped). If `NULL` (default), the day after the
#'   latest event is used. Coerced as for `start_date`. **The modelled study
#'   period always ends the day after the last event kept** (`TT`), whatever
#'   `end_date` is. A later `end_date` is used only to label outputs, such as
#'   the printed study period and the date axis of [plot_event_summary()],
#'   and a message says so.
#' @param crs Target projected CRS (anything accepted by [sf::st_crs()], e.g.
#'   an EPSG code) that both `points` and `region` are transformed to before
#'   any distance/kernel computation. If `NULL` (default), the CRS of
#'   `region` is used, and `points` is transformed to match. Both the
#'   resulting `points` and `region` CRS must be projected (linear units),
#'   not geographic (long/lat) — an informative error is raised otherwise,
#'   and a warning if the projected units are something other than metres
#'   (see the Units section below).
#' @param time_mark_step Step (in days) used to build the `time.marks` grid
#'   returned in the output; default `0.05` matches the original script.
#'
#' @section Units:
#' Every distance in this package — kernel bandwidths, `cell_size`, trigger
#' window distances, triggering-kernel grid extents — is expressed in the
#' linear units of the working CRS, and every default was chosen for a
#' **metre-based** projection. A geographic (longitude/latitude) CRS is
#' rejected outright, because degree-valued distances would be silently
#' meaningless. A projected CRS in feet or kilometres is allowed but warned
#' about, since the defaults will then be wrong by a large factor: usable, but
#' only if you set every distance argument explicitly.
#'
#' @return An object of class `"semihawkes_events"`, a list with elements:
#'   \describe{
#'     \item{events}{A `data.frame` with columns `coorx`, `coory`, `date`,
#'       `days` (numeric days since `start_date`), and `month` (calendar
#'       month, `0`-indexed, i.e. January = `0`).}
#'     \item{region}{`region`, transformed to the working CRS.}
#'     \item{crs}{The working `sf::st_crs()` object.}
#'     \item{start_date, end_date}{The resolved study-period bounds.
#'       `end_date` labels outputs; the model uses `TT`.}
#'     \item{TT}{The modelled study horizon in days, `max(days) + 1`: the
#'       period runs from `start_date` to the day after the last event.}
#'     \item{time.marks}{Numeric grid `seq(0, TT, time_mark_step)`.}
#'     \item{n_input, n_outside_region, n_outside_daterange, n_kept}{Counters
#'       describing how many input rows were dropped and why, for a quick
#'       sanity check on any new dataset.}
#'   }
#'
#' @examples
#' # `sim_events` / `sim_region` are simulated data shipped with the
#' # package; see `?sim_events`.
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' print(ev)
#'
#' # Restricting the study period drops events outside it, and reports how
#' # many were dropped and why.
#' ev_2022 <- prepare_events(sim_events, sim_region, date_col = "date",
#'                           start_date = "2022-01-01",
#'                           end_date   = "2023-01-01")
#' print(ev_2022)
#' @export
prepare_events <- function(points, region, date_col,
                           start_date = NULL, end_date = NULL,
                           crs = NULL, time_mark_step = 0.05) {

  .check_sf_points(points, "points")
  .check_sf_polygon(region, "region")

  if (!is.character(date_col) || length(date_col) != 1L) {
    stop("`date_col` must be a single column name (character scalar).",
         call. = FALSE)
  }
  if (!date_col %in% names(points)) {
    stop(sprintf("`date_col` ('%s') is not a column of `points`. Available columns: %s.",
                 date_col, paste(names(points), collapse = ", ")), call. = FALSE)
  }
  n_input <- nrow(points)

  ## --- 1. Harmonize CRS -----------------------------------------------
  # Only reproject what is not already in the target CRS. Transforming to the
  # same CRS is a no-op for most systems, but GDAL/PROJ cannot build any
  # transformation for a local engineering CRS (e.g. the "Local Cartesian
  # coordinates" of the bundled `crimes` data), not even to itself, and
  # `sf::st_transform()` then fails with "transformation not available".
  target_crs <- if (is.null(crs)) sf::st_crs(region) else sf::st_crs(crs)
  if (!isTRUE(sf::st_crs(points) == target_crs)) {
    points <- sf::st_transform(points, crs = target_crs)
  }
  if (!isTRUE(sf::st_crs(region) == target_crs)) {
    region <- sf::st_transform(region, crs = target_crs)
  }

  .check_projected_crs(points, "points (after transforming to `crs`)")
  .check_projected_crs(region, "region (after transforming to `crs`)")
  .check_crs_units(region, "region")

  ## --- 2. Parse dates and apply the date-range filter ------------------
  raw_dates <- points[[date_col]]
  event_dates <- if (inherits(raw_dates, "Date")) {
    raw_dates
  } else if (inherits(raw_dates, "POSIXt")) {
    .local_date(raw_dates)
  } else {
    as.Date(as.character(raw_dates))
  }
  if (anyNA(event_dates)) {
    warning(sprintf("%d event date(s) could not be parsed and will be dropped.",
                    sum(is.na(event_dates))), call. = FALSE)
  }

  end_date_given <- !is.null(end_date)
  start_date <- if (is.null(start_date)) min(event_dates, na.rm = TRUE) else .local_date(start_date)
  end_date   <- if (is.null(end_date)) max(event_dates, na.rm = TRUE) + 1 else .local_date(end_date)
  if (end_date <= start_date) {
    stop("`end_date` must be strictly after `start_date`.", call. = FALSE)
  }

  in_daterange <- !is.na(event_dates) & event_dates >= start_date & event_dates < end_date
  n_outside_daterange <- sum(!in_daterange)

  points <- points[in_daterange, ]
  event_dates <- event_dates[in_daterange]

  if (nrow(points) == 0) {
    stop("No events remain after filtering by `region` and the date range. ",
         "Check that `points` and `region` share the same area/CRS and that ",
         "`start_date`/`end_date` overlap the data.", call. = FALSE)
  }

  ## --- 3. Build the flat coordinate / time-index structure --------------
  coords <- sf::st_coordinates(points)
  events <- data.frame(
    coorx = coords[, "X"],
    coory = coords[, "Y"],
    date  = event_dates,
    days  = as.numeric(event_dates - start_date),
    month = as.numeric(format(event_dates, "%m")) - 1
  )

  ## --- 4. Keep only events inside `region` ------------------------------
  w <- spatstat.geom::as.owin(region)
  inside <- spatstat.geom::inside.owin(x = events$coorx, y = events$coory, w = w)
  n_outside_region <- sum(!inside)
  events <- events[inside, ]

  if (nrow(events) == 0) {
    stop("No events remain after filtering by `region` and the date range. ",
         "Check that `points` and `region` share the same area/CRS and that ",
         "`start_date`/`end_date` overlap the data.", call. = FALSE)
  }

  ## --- 5. Order by time and build derived quantities ---------------------
  events <- events[order(events$days), ]
  rownames(events) <- NULL

  TT <- range(events$days)[2] + 1
  time.marks <- seq(0, TT, time_mark_step)

  # The modelled study period ends the day after the last event kept. A
  # default `end_date` follows it (events outside `region` may have been
  # dropped after the date filter); a later `end_date` given by the user is
  # kept, but only to label outputs.
  model_end <- start_date + TT
  if (!end_date_given) {
    end_date <- model_end
  } else if (end_date > model_end) {
    message(sprintf(
      paste0("The modelled study period ends on %s, the day after the last ",
             "event. `end_date` (%s) is used only to select events and to ",
             "label outputs."),
      format(model_end), format(end_date)))
  }

  out <- list(
    events = events,
    region = region,
    crs = target_crs,
    start_date = start_date,
    end_date = end_date,
    TT = TT,
    time.marks = time.marks,
    n_input = n_input,
    n_outside_region = n_outside_region,
    n_outside_daterange = n_outside_daterange,
    n_kept = nrow(events)
  )
  class(out) <- "semihawkes_events"
  out
}

#' @export
print.semihawkes_events <- function(x, ...) {
  cat("<semihawkes_events>\n")
  cat(sprintf("  Study period : %s to %s (%sTT = %.1f days)\n",
              x$start_date, x$end_date, .modelled_to(x$start_date, x$end_date, x$TT),
              x$TT))
  cat(sprintf("  Events kept  : %d / %d input rows\n", x$n_kept, x$n_input))
  cat(sprintf("  Dropped      : %d outside region, %d outside date range\n",
              x$n_outside_region, x$n_outside_daterange))
  cat(sprintf("  CRS          : %s\n", .crs_label(x$crs)))
  invisible(x)
}

# "modelled to <date>, " when the modelled period (which ends the day after
# the last event) stops before the labelled `end_date`; "" otherwise.
.modelled_to <- function(start_date, end_date, TT) {
  model_end <- start_date + TT
  if (end_date > model_end) sprintf("modelled to %s, ", format(model_end)) else ""
}

#' @export
summary.semihawkes_events <- function(object, ...) {
  print(object)
  cat("\n")
  summary(object$events)
}

# Calendar date of a Date, date-time or date string. A date-time is converted
# in its own time zone: its `tzone` attribute, or the session's time zone if
# it has none. as.Date() on a POSIXct uses UTC instead, which moves evening
# events west of Greenwich, and early-morning events east of it, to the
# neighbouring day.
.local_date <- function(x) {
  if (inherits(x, "POSIXt")) as.Date(format(x, "%Y-%m-%d")) else as.Date(x)
}

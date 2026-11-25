#' Day-of-week and calendar-month marks for a set of day offsets
#'
#' Converts numeric day offsets (days since `start_date`) into the two
#' periodic marks used by [temporal_base_weekly()] and
#' [temporal_base_monthly()]: day-of-week (`0` = Sunday, matching
#' `as.POSIXlt()$wday`) and calendar month (`0`-indexed, matching the
#' `month` column built in `prepare_events()`). Needed because `weekly$fun`
#' and `monthly$fun` are periodic in *calendar* time, not in days-since-start,
#' so callers must not pass raw `days` offsets to them directly.
#'
#' @param days Numeric vector of day offsets since `start_date`.
#' @param start_date The study period's `start_date` (a `Date`).
#' @return A list with `wday` and `month`, both integer vectors the same
#'   length as `days`.
#' @keywords internal
#' @noRd
.calendar_marks <- function(days, start_date) {
  dates <- start_date + days
  lt <- as.POSIXlt(dates)
  list(
    wday  = lt$wday,
    month = lt$mon
  )
}

# Non-parametric temporal background: weekly, monthly (calendar-month) and
# long-term trend base values, generalizing the "BASE VALUES" section of the
# original script (both the Loop-0 initial guess and the Loop-k E-step
# re-estimation, which reuse the same `hist.weighted()` + `ker.smooth.fft()`
# / weighted-Gaussian machinery with different weights).
#
#' Count how often each periodic "mark" occurs over the full study calendar
#'
#' Generalizes the day-of-week / day-of-month count-correction weights from
#' the original script: instead of an algebraic approximation, this walks
#' every calendar day in `[start_date, end_date)` and tallies how often each
#' value of `mark_fun(date)` occurs, so events falling on an
#' over-/under-represented weekday or month (because the study period isn't
#' an exact multiple of the period) get corrected weights.
#'
#' @param start_date,end_date Study-period bounds (`end_date` exclusive).
#' @param mark_fun A function mapping a vector of `Date`s to periodic mark
#'   values (e.g. `function(d) as.POSIXlt(d)$wday` for day-of-week).
#'
#' @return A named numeric vector: counts per distinct mark value, named by
#'   the (character) mark value.
#' @keywords internal
#' @noRd
.calendar_mark_counts <- function(start_date, end_date, mark_fun) {
  all_dates <- seq(as.Date(start_date), as.Date(end_date) - 1, by = "day")
  table(mark_fun(all_dates))
}

#' Per-event count-correction weight for a periodic mark
#'
#' `1 / count` of each event's mark, using `.calendar_mark_counts()` so that
#' marks that occur more often within `[start_date, end_date)` (e.g. because
#' the horizon isn't an exact number of weeks/years) don't bias the
#' non-parametric density estimate.
#'
#' @param start_date,end_date Study-period bounds (`end_date` exclusive).
#' @param mark_fun A function mapping a vector of `Date`s to periodic mark
#'   values.
#' @param dates `Date` vector of event dates to weight (a subset of, or
#'   overlapping with, `[start_date, end_date)`).
#' @return A numeric vector, the same length as `dates`.
#' @keywords internal
#' @noRd
.calendar_mark_weights <- function(dates, start_date, end_date, mark_fun) {
  counts <- .calendar_mark_counts(start_date, end_date, mark_fun)
  marks <- as.character(mark_fun(dates))
  as.numeric(1 / counts[marks])
}

#' Build a periodic evaluator function from a base-value curve
#'
#' Wraps [stats::approxfun()] so it can be evaluated at any `x` by first
#' reducing `x` modulo `period` (matching the original script's
#' `x - as.integer(x / period) * period`, which is exact for the
#' non-negative day offsets used throughout this package).
#'
#' @param breaks,basevalue Grid and values defining the periodic curve.
#' @param period Numeric period (`7` for weekly, `12` for monthly).
#' @return A vectorized function of `x`.
#' @keywords internal
#' @noRd
.periodic_function <- function(breaks, basevalue, period) {
  f <- stats::approxfun(breaks, basevalue, yleft = 0, yright = 0)
  function(x) f(x - as.integer(x / period) * period)
}

#' Weekly non-parametric temporal base value
#'
#' Estimates a periodic (7-day) background multiplier from event
#' days-of-week, following the weekly block of the original
#' Zhuang & Mateu (2019) style script: bin events by day-of-week with
#' [hist.weighted()], smooth with [ker.smooth.fft()], and normalize to a
#' mean of 1. Generalizes both the Loop-0 initial guess (`resp_weight = NULL`)
#' and the Loop-k E-step re-estimation (`resp_weight` = the EM responsibility
#' weight for the weekly component at each event, e.g. `weekly.fun(days) *
#' mub.events / lambda.at.events` in the original script) into one function.
#'
#' @param events A `"semihawkes_events"` object from [prepare_events()].
#' @param resp_weight Optional numeric vector, one entry per event in
#'   `events$events`, e.g. an EM E-step responsibility weight. If `NULL`
#'   (default), only the day-count correction weight is used (the Loop-0
#'   initial-guess behaviour).
#' @param breaks Numeric grid of day-of-week bin edges; default
#'   `seq(0, 7, 0.05)` (`0` = Sunday), matching the original script.
#' @param bandwidth Gaussian smoothing bandwidth passed to
#'   [ker.smooth.fft()]; default `1`.
#'
#' @return A list of class `"semihawkes_temporal_base"` with elements
#'   `breaks`, `basevalue` (normalized to mean 1), `period` (`7`), and `fun`
#'   (a vectorized function of `days` returning the periodic multiplier).
#' @references Zhuang & Mateu (2019); see [spHawkes-package].
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' weekly <- temporal_base_weekly(ev)
#' weekly
#' # Multiplier for each day of the week (0 = Sunday)
#' round(weekly$fun(0:6), 2)
#' @export
temporal_base_weekly <- function(events, resp_weight = NULL,
                                 breaks = seq(0, 7, 0.05), bandwidth = 1) {
  stopifnot(inherits(events, "semihawkes_events"))
  dates <- events$events$date
  wday <- as.POSIXlt(dates)$wday # 0 = Sunday, ..., 6 = Saturday

  count_weight <- .calendar_mark_weights(
    dates, events$start_date, events$end_date,
    mark_fun = function(d) as.POSIXlt(d)$wday
  )
  total_weight <- count_weight * (if (is.null(resp_weight)) 1 else resp_weight)

  hw <- hist.weighted(wday, total_weight, breaks = breaks)
  basevalue <- ker.smooth.fft(hw$mids, hw$density, bandwidth) # length == length(breaks)
  basevalue <- basevalue / mean(basevalue)

  out <- list(breaks = breaks, basevalue = basevalue, period = 7,
              fun = .periodic_function(breaks, basevalue, period = 7))
  class(out) <- "semihawkes_temporal_base"
  out
}

#' Monthly non-parametric temporal base value
#'
#' Estimates a periodic (12-month) background multiplier from each event's
#' calendar month (`events$events$month`, `0` = January), the same way as
#' [temporal_base_weekly()] but with a 12-month period. See that function's
#' documentation for the `resp_weight` argument and the Loop-0 vs. Loop-k
#' generalization.
#'
#' @inheritParams temporal_base_weekly
#' @param breaks Numeric grid of month-of-year bin edges; default
#'   `seq(0, 12, 0.05)`.
#'
#' @return A list of class `"semihawkes_temporal_base"`; see
#'   [temporal_base_weekly()].
#' @references Zhuang & Mateu (2019); see [spHawkes-package].
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' monthly <- temporal_base_monthly(ev)
#' # Multiplier for each calendar month (0 = January)
#' round(monthly$fun(0:11), 2)
#' @export
temporal_base_monthly <- function(events, resp_weight = NULL,
                                  breaks = seq(0, 12, 0.05), bandwidth = 1) {
  stopifnot(inherits(events, "semihawkes_events"))
  dates <- events$events$date
  month_mark <- events$events$month # 0-indexed calendar month, from prepare_events()

  count_weight <- .calendar_mark_weights(
    dates, events$start_date, events$end_date,
    mark_fun = function(d) as.numeric(format(d, "%m")) - 1
  )
  total_weight <- count_weight * (if (is.null(resp_weight)) 1 else resp_weight)

  hw <- hist.weighted(month_mark, total_weight, breaks = breaks)
  basevalue <- ker.smooth.fft(hw$mids, hw$density, bandwidth) # length == length(breaks)
  basevalue <- basevalue / mean(basevalue)

  out <- list(breaks = breaks, basevalue = basevalue, period = 12,
              fun = .periodic_function(breaks, basevalue, period = 12))
  class(out) <- "semihawkes_temporal_base"
  out
}

#' Long-term trend temporal base value
#'
#' Estimates a non-periodic, long-term background multiplier over the whole
#' study horizon. With `resp_weight = NULL` (default) this returns a flat
#' trend (`basevalue = 1` everywhere) - the Loop-0 initial guess in the
#' original script, where the weighted-smoothing loop was disabled and the
#' trend was left flat until the first E-step. With `resp_weight` supplied
#' (one EM responsibility weight per event, e.g. `trend.fun(days) *
#' mub.events / lambda.at.events`), it reproduces the Loop-k weighted Gaussian
#' smoother with boundary (edge-effect) correction via [stats::pnorm()].
#'
#' Events fall on whole days, so the weights of events on the same day are
#' added up first, and the smoothed curve is evaluated on a grid of one day,
#' or `bandwidth / 100` days if that is finer, then interpolated linearly to
#' `events$time.marks`. That changes the curve by at most about
#' `(step / bandwidth)^2 / 8` of its value: 6 parts in a million for the
#' default bandwidth. When the grid would be as fine as the time marks, the
#' curve is evaluated at the time marks directly.
#'
#' @inheritParams temporal_base_weekly
#' @param bandwidth Gaussian smoothing bandwidth (in days); default `150`,
#'   matching the original script.
#' @param chunk_size Number of distinct event days processed per chunk when
#'   `resp_weight` is supplied, to bound peak memory use for long study
#'   periods; default `500`.
#'
#' @return A list of class `"semihawkes_temporal_base"` with elements
#'   `breaks` (`events$time.marks`), `basevalue` (normalized to mean 1),
#'   `period` (`NA`, i.e. non-periodic), and `fun` (a vectorized function of
#'   `days`, `0` outside `[0, TT]`).
#' @references Zhuang & Mateu (2019); see [spHawkes-package].
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date",
#'                      end_date = "2022-01-01")
#' # Without weights: the flat initial guess used before the first EM step
#' flat <- temporal_base_trend(ev)
#' range(flat$basevalue)
#' # With weights (here every event weighted equally): a smoothed trend
#' trend <- temporal_base_trend(ev, resp_weight = rep(1, nrow(ev$events)))
#' plot(trend$breaks, trend$basevalue, type = "l",
#'      xlab = "Days since start", ylab = "Trend multiplier")
#' @export
temporal_base_trend <- function(events, resp_weight = NULL, bandwidth = 150,
                                chunk_size = 500) {
  stopifnot(inherits(events, "semihawkes_events"))
  time_marks <- events$time.marks

  if (is.null(resp_weight)) {
    basevalue <- rep(1, length(time_marks))
  } else {
    event_days <- events$events$days
    n <- length(event_days)
    stopifnot(length(resp_weight) == n)

    # Each event's weight, divided by its kernel's mass inside [0, TT]
    # (the edge correction), then added up over events on the same day.
    denom <- stats::pnorm(events$TT, event_days, bandwidth) -
      stats::pnorm(0, event_days, bandwidth)
    days <- sort(unique(event_days))
    day_weight <- as.numeric(rowsum(resp_weight / denom, match(event_days, days)))

    # Evaluate on a grid no finer than the curve needs (see Details).
    step <- min(1, bandwidth / 100)
    mark_step <- if (length(time_marks) > 1) time_marks[2] - time_marks[1] else Inf
    grid <- if (step <= mark_step) {
      time_marks
    } else {
      seq(0, max(time_marks) + step, by = step)
    }

    on_grid <- numeric(length(grid))
    for (ix in split(seq_along(days), ceiling(seq_along(days) / chunk_size))) {
      dens <- stats::dnorm(outer(days[ix], grid, FUN = "-"), 0, bandwidth)
      on_grid <- on_grid + colSums(day_weight[ix] * dens)
    }
    basevalue <- if (identical(grid, time_marks)) {
      on_grid
    } else {
      stats::approx(grid, on_grid, xout = time_marks)$y
    }
    basevalue <- basevalue / mean(basevalue)
  }

  out <- list(breaks = time_marks, basevalue = basevalue, period = NA_real_,
              fun = stats::approxfun(time_marks, basevalue, yleft = 0, yright = 0))
  class(out) <- "semihawkes_temporal_base"
  out
}

#' @export
print.semihawkes_temporal_base <- function(x, ...) {
  cat("<semihawkes_temporal_base>\n")
  cat(sprintf("  Grid points : %d (range %.3g to %.3g)\n",
              length(x$breaks), min(x$breaks), max(x$breaks)))
  cat(sprintf("  Period      : %s\n", if (is.na(x$period)) "none (trend)" else x$period))
  cat(sprintf("  Base value  : mean %.3g, range [%.3g, %.3g]\n",
              mean(x$basevalue), min(x$basevalue), max(x$basevalue)))
  invisible(x)
}

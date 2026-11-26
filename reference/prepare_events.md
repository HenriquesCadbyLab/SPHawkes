# Prepare a spatio-temporal event dataset for semi-parametric Hawkes fitting

Validates, cleans and reshapes a set of georeferenced events (e.g. crime
records) and a study-region polygon into the flat structure used by the
rest of the `spHawkes` EM pipeline: projected `x`/`y` coordinates, a
numeric time index in days since `start_date`, and events restricted to
`region` and to `[start_date, end_date)`.

## Usage

``` r
prepare_events(
  points,
  region,
  date_col,
  start_date = NULL,
  end_date = NULL,
  crs = NULL,
  time_mark_step = 0.05
)
```

## Arguments

- points:

  An `sf` object with `POINT` geometry (one row per event).

- region:

  An `sf`/`sfc` object with `POLYGON`/`MULTIPOLYGON` geometry describing
  the study area (e.g. a city or district boundary).

- date_col:

  Name of the column in `points` holding the event date (character
  scalar). The column may be `Date`, `POSIXct`/`POSIXlt`, or a string
  parseable by [`as.Date()`](https://rdrr.io/r/base/as.Date.html).
  Date-times are converted to calendar dates in their own time zone
  (their `tzone` attribute, or the session's time zone if they have
  none), so an event at 20:00 on a Friday in Chicago stays on that
  Friday.

- start_date:

  Start of the study period. If `NULL` (default), the earliest event
  date is used. Coerced with
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html), or, for a
  date-time, converted in its own time zone as for `date_col`.

- end_date:

  Last date for selecting events, exclusive (events with
  `date >= end_date` are dropped). If `NULL` (default), the day after
  the latest event is used. Coerced as for `start_date`. **The modelled
  study period always ends the day after the last event kept** (`TT`),
  whatever `end_date` is. A later `end_date` is used only to label
  outputs, such as the printed study period and the date axis of
  [`plot_event_summary()`](https://henriquescadbylab.github.io/spHawkes/reference/plot_event_summary.md),
  and a message says so.

- crs:

  Target projected CRS (anything accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html),
  e.g. an EPSG code) that both `points` and `region` are transformed to
  before any distance/kernel computation. If `NULL` (default), the CRS
  of `region` is used, and `points` is transformed to match. Both the
  resulting `points` and `region` CRS must be projected (linear units),
  not geographic (long/lat) — an informative error is raised otherwise,
  and a warning if the projected units are something other than metres
  (see the Units section below).

- time_mark_step:

  Step (in days) used to build the `time.marks` grid returned in the
  output; default `0.05` matches the original script.

## Value

An object of class `"semihawkes_events"`, a list with elements:

- events:

  A `data.frame` with columns `coorx`, `coory`, `date`, `days` (numeric
  days since `start_date`), and `month` (calendar month, `0`-indexed,
  i.e. January = `0`).

- region:

  `region`, transformed to the working CRS.

- crs:

  The working
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  object.

- start_date, end_date:

  The resolved study-period bounds. `end_date` labels outputs; the model
  uses `TT`.

- TT:

  The modelled study horizon in days, `max(days) + 1`: the period runs
  from `start_date` to the day after the last event.

- time.marks:

  Numeric grid `seq(0, TT, time_mark_step)`.

- n_input, n_outside_region, n_outside_daterange, n_kept:

  Counters describing how many input rows were dropped and why, for a
  quick sanity check on any new dataset.

## Units

Every distance in this package — kernel bandwidths, `cell_size`, trigger
window distances, triggering-kernel grid extents — is expressed in the
linear units of the working CRS, and every default was chosen for a
**metre-based** projection. A geographic (longitude/latitude) CRS is
rejected outright, because degree-valued distances would be silently
meaningless. A projected CRS in feet or kilometres is allowed but warned
about, since the defaults will then be wrong by a large factor: usable,
but only if you set every distance argument explicitly.

## Examples

``` r
# `sim_events` / `sim_region` are simulated data shipped with the
# package; see `?sim_events`.
ev <- prepare_events(sim_events, sim_region, date_col = "date")
print(ev)
#> <semihawkes_events>
#>   Study period : 2021-01-03 to 2024-01-01 (TT = 1093.0 days)
#>   Events kept  : 375 / 375 input rows
#>   Dropped      : 0 outside region, 0 outside date range
#>   CRS          : EPSG:32717

# Restricting the study period drops events outside it, and reports how
# many were dropped and why.
ev_2022 <- prepare_events(sim_events, sim_region, date_col = "date",
                          start_date = "2022-01-01",
                          end_date   = "2023-01-01")
print(ev_2022)
#> <semihawkes_events>
#>   Study period : 2022-01-01 to 2023-01-01 (TT = 365.0 days)
#>   Events kept  : 113 / 375 input rows
#>   Dropped      : 0 outside region, 262 outside date range
#>   CRS          : EPSG:32717
```

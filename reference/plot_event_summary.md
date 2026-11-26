# Plot a summary of a spatio-temporal event data set

Creates a three-panel diagnostic display containing the event locations,
gridded spatial counts, and cumulative event count.

## Usage

``` r
plot_event_summary(
  x,
  region = x$region,
  start_date = x$start_date,
  end_date = x$end_date,
  grid_n = 10L
)
```

## Arguments

- x:

  A `semihawkes_events` object, as returned by
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md).

- region:

  An `sf` polygon defining the study region. Defaults to `x$region`.

- start_date, end_date:

  Study bounds; `end_date` is exclusive. They default to `x$start_date`
  and `x$end_date`.

- grid_n:

  Approximate number of square cells across the longest side of the
  region's bounding box.

## Value

A `patchwork` object.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
plot_event_summary(ev)
```

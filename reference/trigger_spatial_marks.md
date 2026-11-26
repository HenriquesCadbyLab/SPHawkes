# Cache per-event spatial edge-correction marks

For each event `i`, shifts the triggering kernel's spatial grid
(`spatial_trigger$x`/`spatial_trigger$y`) by that event's coordinates
and flags which shifted grid cells fall inside the study region, caching
the result to disk (one `.rds` file per event).
[`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md)
reads these back to compute each event's spatial edge-correction factor
and the pooled spatial "repetance" surface.

## Usage

``` r
trigger_spatial_marks(
  events,
  spatial_trigger,
  cache_dir,
  force = FALSE,
  nprocs = 1
)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)).

- spatial_trigger:

  A `semihawkes_trigger_spatial` object (from
  [`trigger_spatial_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_init.md)
  or
  [`trigger_reestimate_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_spatial.md))
  – its `x`/`y` grid defines the marks that get cached.

- cache_dir:

  Directory to write/read `.rds` mark files; created if it doesn't
  exist. Required, no default.

- force:

  If `TRUE`, recompute and overwrite marks that are already cached;
  default `FALSE` (skip events with an existing file, matching the
  original script's `if (!file.exists(fn))` guard).

- nprocs:

  Parallel workers, as for
  [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md):
  a whole number (the default, `1`, runs serially), `"auto"`, or a
  cluster from
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
  which is used and left running.

## Value

`cache_dir`, invisibly.

## Details

A grid cell counts as inside if it intersects the region, as
[`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
decides, so cells on the boundary count as inside. To keep this fast for
many events, the shifted grids of up to 200 events are tested at once,
mostly with spatstat's point-in-polygon test; cells within a hair's
breadth of the boundary (a thousandth of the grid step) are tested with
[`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
itself, so the marks are the same as testing every cell with sf.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
ts <- trigger_spatial_init(x_range = 2000, x_step = 250)
cache_dir <- tempfile("spHawkes_marks_")
trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
length(list.files(cache_dir))  # one file per event
#> [1] 98
```

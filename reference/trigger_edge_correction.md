# Per-event temporal/spatial edge-correction factors for trigger re-estimation

Computes, for every event, how much of its triggering kernels' mass
would have been observable given the study period's end date and the
study region's boundary – the correction factors the original script
calls `excite.temporal.edge.correction` /
`excite.spatial.edge.correction` – plus the pooled "repetance" surfaces
(`temporal.repetance` / `spatial.repetance` in the original) used to
correct the re-estimated kernels for how many events *could* have
produced an observation at each offset, given the same boundaries.

## Usage

``` r
trigger_edge_correction(
  events,
  temporal_trigger,
  spatial_trigger,
  cache_dir,
  marks = NULL
)
```

## Arguments

- events:

  A `semihawkes_events` object.

- temporal_trigger:

  A `semihawkes_trigger_temporal` object.

- spatial_trigger:

  A `semihawkes_trigger_spatial` object.

- cache_dir:

  Directory previously populated by
  [`trigger_spatial_marks()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_marks.md)
  for the same `events`/`spatial_trigger`.

- marks:

  Optional: the marks in `cache_dir`, already read into a list with one
  logical matrix per event (as
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  does once per fit, so they are not read from disk at every iteration).
  `NULL` (the default) reads them from `cache_dir`.

## Value

A list of class `"semihawkes_trigger_edge_correction"` with
`temporal_edge_correction`, `spatial_edge_correction` (numeric vectors,
one value per event), `temporal_repetance` (numeric vector on
`temporal_trigger$base`), `spatial_repetance` (matrix on
`spatial_trigger$x`/`y`), and `temporal_repetance_fun`/
`spatial_repetance_fun` (ready-to-use evaluators).

## Details

Requires
[`trigger_spatial_marks()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_marks.md)
to have already populated `cache_dir` for every event in `events`.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
tt <- trigger_temporal_init()
ts <- trigger_spatial_init(x_range = 2000, x_step = 250)
cache_dir <- tempfile("spHawkes_marks_")
trigger_spatial_marks(ev, ts, cache_dir = cache_dir)
ec <- trigger_edge_correction(ev, tt, ts, cache_dir)
ec
#> <semihawkes_trigger_edge_correction>
#>   98 events
```

# Evaluate background and triggering intensity on a spatial grid over time

Computes, for each day in the study period (or a user-supplied set of
days), the background intensity \\\mu(s,t)\\ and the triggering
intensity \\A \cdot g(t) \cdot f(s)\\ on a regular spatial grid covering
`city`. This is the building block behind exceedance-set and persistence
analyses: it turns a fitted `semihawkes_fit` object into daily snapshots
of both intensity components that can later be aggregated over periods
(yearly, semestral, quarterly, ...) via
[`integrate_intensity_period()`](https://henriquescadbylab.github.io/spHawkes/reference/integrate_intensity_period.md).

## Usage

``` r
background_trigger_grid(
  fit,
  city,
  cell_size = 150,
  days = NULL,
  t.lag = 90,
  s.lag = NULL,
  nprocs = 1
)
```

## Arguments

- fit:

  A `semihawkes_fit` object (or equivalently-shaped list) with
  `$events`, `$mu`, `$A`, `$trend`, `$weekly`, `$monthly`,
  `$background`, `$temporal_trigger`, `$spatial_trigger`.

- city:

  An `sf` polygon giving the study region. Used both to build the
  evaluation grid (from its bounding box) and to mask grid cells outside
  the region (via
  [`spatstat.geom::inside.owin()`](https://rdrr.io/pkg/spatstat.geom/man/inside.owin.html)).

- cell_size:

  Grid spacing, in the same spatial units as `city` (e.g. metres for UTM
  data). Default `150`.

- days:

  Integer vector of days-since-`start_date` at which to evaluate the
  intensity. Defaults to `0:TT`, where
  `TT <- max(fit$events$events$days) + 1`.

- t.lag:

  Only triggering events within this many days before `t` contribute
  (events further in the past are dropped for speed, on the assumption
  their contribution has decayed to ~0). Default `90`.

- s.lag:

  Only triggering events within this spatial distance of the grid's
  bounding box contribute. The default, `NULL`, uses the reach of the
  fitted spatial triggering kernel, beyond which its contribution is
  zero.

- nprocs:

  Number of parallel workers, via `foreach`/`doParallel`. `nprocs = 1`
  (default) runs sequentially with no cluster overhead. With more, the
  workers are stopped on exit and the sequential `foreach` backend is
  registered in their place.

## Value

A list with:

- `days`:

  the days at which intensity was evaluated (matches the `days`
  argument).

- `x`, `y`:

  grid coordinate vectors.

- `inside`:

  logical vector, `TRUE` for grid cells inside `city`.

- `snapshots`:

  a list (one element per day) of lists with `$bg` and `$tr`, each a
  `length(x)` by `length(y)` matrix of intensity values. Cells outside
  `city` are `0`.

## Details

Unlike
[`background_intensity()`](https://henriquescadbylab.github.io/spHawkes/reference/background_intensity.md)
and
[`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md),
which evaluate the model only at observed event locations (for
likelihood/EM purposes), this function evaluates it on an arbitrary
regular grid, independent of the grid resolution
[`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md)
used internally during fitting. You choose the grid resolution via
`cell_size`.

`fit$background$fun`, `fit$trend$fun`, `fit$weekly$fun`,
`fit$monthly$fun`, `fit$temporal_trigger$fun` and
`fit$spatial_trigger$fun` are assumed to be the final, converged fit –
i.e. `fit` is the result of a completed
[`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
run (or the last
[`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md)
in a manual loop), not an initial Loop-0 guess.

## Examples

``` r
# \donttest{
# A quick fit to one year of the bundled simulated data: a coarse grid, a
# coarse triggering-kernel grid and two EM iterations keep it to a few
# seconds. Real analyses need finer settings and more iterations.
cache_dir <- tempfile("spHawkes_marks_")
fit <- fit_semihawkes(
  sim_events, sim_region,
  date_col  = "date",
  end_date  = "2022-01-01",
  cell_size = 1000,
  cache_dir = cache_dir,
  trigger_spatial_args = list(x_range = 2000, x_step = 250),
  trigger_window_dist  = 2000,
  n_iter  = 2,
  verbose = FALSE,
  seed    = 1
)
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
# Evaluate on a handful of days rather than the whole study window
grid <- background_trigger_grid(fit, sim_region, cell_size = 500,
                                days = c(0, 100, 200))
length(grid$snapshots)
#> [1] 3
dim(grid$snapshots[[1]]$bg)
#> [1] 17 17
# }
```

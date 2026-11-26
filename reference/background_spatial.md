# Non-parametric spatial background surface

Estimates the spatial component of the background intensity on a regular
grid covering the bounding box of `events$region`, by summing each
event's isotropic Gaussian kernel (using its adaptive `bandwidth` and
edge-correction `bg.integral`, see
[`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md)
and
[`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md))
evaluated at every grid cell, then mean-normalizing over the cells that
fall inside `region`.

## Usage

``` r
background_spatial(events, cell_size, resp_weight = NULL, chunk_size = 200)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)),
  already passed through
  [`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md)
  and
  [`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md)
  (i.e. `events$events$bandwidth` and `events$events$bg.integral` must
  exist).

- cell_size:

  Grid cell size, in the same linear units as the working CRS (e.g.
  metres for a UTM projection). There is deliberately no default: pick a
  value appropriate to your CRS units and the spatial scale of your data
  (the original script used `150` for a dataset projected in metres).

- resp_weight:

  Optional numeric vector, one value per row of `events$events`, giving
  the EM E-step background-responsibility weight for each event. If
  `NULL` (the default), every event gets weight `1` (the initial-guess
  behaviour).

- chunk_size:

  Number of events processed per vectorized chunk, to cap memory use for
  large datasets and/or fine grids (each chunk builds a
  `n_grid_cells x chunk_size` matrix internally). Default `200`.

## Value

A list of class `"semihawkes_spatial_base"` with elements:

- x, y:

  The grid's coordinate sequences.

- basevalue:

  `length(x) x length(y)` matrix, mean-normalized over cells with
  `marks > 0`.

- marks:

  `length(x) x length(y)` integer matrix, `1` where the cell centre
  falls inside `events$region`, `0` otherwise.

- cell_size:

  The `cell_size` used.

- fun:

  A function `fun(x, y)` that bilinearly interpolates `basevalue` at
  arbitrary coordinates via
  [`fields::interp.surface()`](https://rdrr.io/pkg/fields/man/interp.surface.html),
  returning `NA` outside the grid.

## Details

As with
[`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md),
[`temporal_base_monthly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_monthly.md)
and
[`temporal_base_trend()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_trend.md),
`resp_weight = NULL` reproduces the EM initial-guess step (every event
contributes with weight 1), and supplying `resp_weight` (the E-step
background-responsibility vector \\\varphi_i\\, one value per event)
reproduces the weighted re-estimation step.

## References

Zhuang, J. and Mateu, J. (2019). A semiparametric spatiotemporal
Hawkes-type point process model with periodic background for crime data.
*Journal of the Royal Statistical Society Series A*, 182(3), 919-942.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
bg <- background_spatial(ev, cell_size = 1000)
bg
#> <semihawkes_spatial_base>
#>   grid: 9 x 9 cells (cell_size = 1000), 81 inside `region`
#>   basevalue: mean-normalized over in-region cells (mean = 1.0000)
```

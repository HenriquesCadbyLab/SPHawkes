# Integrate gridded intensity over one time period

Numerically integrates background and triggering intensity snapshots
returned by background_trigger_grid(), using the trapezoidal rule.

## Usage

``` r
integrate_intensity_period(grid, start, end)
```

## Arguments

- grid:

  A list returned by background_trigger_grid().

- start:

  Initial model time, in days.

- end:

  Final model time, in days. Must be greater than start.

## Value

A list with days, x, y, inside and snapshots, compatible with
plot_intensity_snapshot(). The single snapshot contains temporally
integrated bg and tr matrices. Its time identifier is start. Additional
fields statistic, period and integration describe the aggregation.

## Details

Integration is performed independently at each spatial grid location.
Actual time differences are used, so equally spaced times are not
required.

Both period boundaries must be available in grid\$days. No temporal
interpolation or extrapolation is performed.

This function integrates sampled intensity surfaces. It does not
reevaluate the fitted model or implement event-aware quadrature.
Accuracy should be assessed by refining the evaluation times,
particularly near triggering discontinuities.

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
grid <- background_trigger_grid(fit, sim_region, cell_size = 500,
                                days = 0:30)
january <- integrate_intensity_period(grid, start = 0, end = 30)
str(january, max.level = 1)
#> List of 8
#>  $ days       : num 0
#>  $ x          : num [1:17] 0 500 1000 1500 2000 2500 3000 3500 4000 4500 ...
#>  $ y          : num [1:17] 0 500 1000 1500 2000 2500 3000 3500 4000 4500 ...
#>  $ inside     : logi [1:289] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ snapshots  :List of 1
#>  $ statistic  : chr "time_integral"
#>  $ period     : Named num [1:2] 0 30
#>   ..- attr(*, "names")= chr [1:2] "start" "end"
#>  $ integration:List of 3
# }
```

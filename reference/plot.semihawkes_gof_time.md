# Plot a time-transformed goodness-of-fit diagnostic

Plots the observed cumulative event count against the model-implied
transformed time, together with pointwise KS-style acceptance bands.

## Usage

``` r
# S3 method for class 'semihawkes_gof_time'
plot(x, ...)
```

## Arguments

- x:

  An object of class `semihawkes_gof_time` returned by
  [`gof_time_transform()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_time_transform.md).

- ...:

  Additional arguments passed to ggplot2 functions.

## Value

A `ggplot` object.

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
plot(gof_time_transform(fit))

# }
```

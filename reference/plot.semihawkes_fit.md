# Plot method for a fitted `semihawkes_fit` object

Replaces the script's scattered
[`png()`](https://rdrr.io/r/grDevices/png.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/`ggplot()`/
[`filled.contour()`](https://rdrr.io/r/graphics/filled.contour.html)
calls and the final `patchwork` combination with a single method
covering weekly/monthly/trend, the spatial background, the two
triggering kernels, and the EM convergence history (log-likelihood, mu,
A, lambda).

## Usage

``` r
# S3 method for class 'semihawkes_fit'
plot(x, which = "history", ...)
```

## Arguments

- x:

  A `semihawkes_fit` object.

- which:

  One of `"weekly"`, `"monthly"`, `"trend"`, `"background"`,
  `"trigger_temporal"`, `"trigger_spatial"`, `"history"`.

- ...:

  Currently unused.

## Value

`x`, invisibly. Called for its side effect of drawing the plot.

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
plot(fit, which = "weekly")

plot(fit, which = "trigger_temporal")

plot(fit, which = "background")

# }
```

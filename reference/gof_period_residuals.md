# Spatial residuals by time period

Computes spatial raw residuals separately for each calendar year or
calendar month in the fitted observation window.

## Usage

``` r
gof_period_residuals(fit, period = c("year", "month"))
```

## Arguments

- fit:

  A fitted object of class `semihawkes_fit` returned by
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).

- period:

  Either `"year"` or `"month"`.

## Value

A named list of `semihawkes_gof_spatial` objects, one for each calendar
year (named like `"2021"`) or month (`"2021-03"`) that overlaps the
fitted observation window, including any without events.

## Details

Each window runs from the first to the last day of its year or month,
cut to the fitted observation window at either end. Together the windows
cover the observation window exactly once, so the observed and the
expected counts both add up to their totals over the whole window.

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
res <- gof_period_residuals(fit, period = "month")
plot(res)

# }
```

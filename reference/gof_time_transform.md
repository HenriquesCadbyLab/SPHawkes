# Time-transformed goodness-of-fit diagnostic

Computes the model-implied cumulative expected count over time,
\$\$\tau(t) = \int_0^t \Lambda(u) \\ du,\$\$ where \\\Lambda(t)\\ is the
fitted intensity integrated over the study region, and compares it with
the observed cumulative event count.

## Usage

``` r
gof_time_transform(fit, dt = 1, level = 0.95)
```

## Arguments

- fit:

  A fitted object of class `semihawkes_fit` returned by
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).

- dt:

  Temporal grid spacing, in days, used to approximate the time integral.

- level:

  Coverage level of the pointwise KS-style bands.

## Value

An object of class `semihawkes_gof_time`: a list containing `grid`,
`events`, `bands` and `level`.

## Details

Under a correctly specified temporal intensity, the transformed event
times are approximately a unit-rate Poisson process. The function also
returns pointwise KS-style acceptance bands for the observed cumulative
count.

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
gof <- gof_time_transform(fit, dt = 1)
plot(gof)

# }
```

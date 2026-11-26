# Total model intensity integrated over space

Computes the fitted semi-parametric Hawkes intensity integrated over the
study region at a vector of times, \$\$\Lambda(t) = \int_W \lambda(s, t)
\\ ds.\$\$

## Usage

``` r
intensity_at_times(fit, times, exclude_self = TRUE)
```

## Arguments

- fit:

  A fitted object of class `semihawkes_fit` returned by
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).

- times:

  Numeric vector of times, in days, at which to evaluate the spatially
  integrated intensity. Times must lie within the fitted observation
  window.

- exclude_self:

  Logical; if `TRUE` (default), an event does not contribute to its own
  triggering intensity.

## Value

A numeric vector, of the same length as `times`, containing the model
intensity integrated over the study region at each time.

## Details

The background contribution is \$\$\mu \mu_T(t) \mu_w(t) \mu_m(t) \int_W
\mu_b(s) \\ ds,\$\$ while the triggering contribution is \$\$A
\sum\_{i:t_i \< t} g(t - t_i) \int_W h(s - s_i) \\ ds.\$\$

Spatial integrals use the same background grid and edge-correction marks
as
[`background_intensity()`](https://henriquescadbylab.github.io/spHawkes/reference/background_intensity.md)
and
[`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md),
so the results are numerically consistent with `fit$lambda_at_all`.

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
# Expected events per day across the study region, on three days
intensity_at_times(fit, times = c(30, 180, 330))
#> [1] 0.1843796 0.2940002 0.1952269
# }
```

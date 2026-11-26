# Single EM step: E-step responsibilities + M-step update of A, mu, and all non-parametric background/triggering components

Single EM step: E-step responsibilities + M-step update of A, mu, and
all non-parametric background/triggering components

## Usage

``` r
em_step(
  events,
  background,
  trend,
  weekly,
  monthly,
  temporal_trigger,
  spatial_trigger,
  mu,
  A,
  pairs,
  cache_dir,
  cell_size,
  nprocs = 1,
  exclude_self = TRUE,
  bgint = NULL,
  trig = NULL,
  marks = NULL
)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md),
  already passed through
  [`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md)
  /
  [`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md)).

- background:

  A `semihawkes_spatial_base` object from the *previous* iteration's
  [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md)
  call (only its grid/marks are reused as the spatial domain; a fresh
  weighted surface is computed inside this step).

- trend, weekly, monthly:

  The *previous* iteration's `semihawkes_temporal_base` objects
  (weekly/monthly/trend).

- temporal_trigger, spatial_trigger:

  The *previous* iteration's triggering-kernel objects.

- mu, A:

  The *previous* iteration's scalar parameter estimates.

- pairs:

  A
  [`trigger_pairs()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pairs.md)
  index (built once, reused every iteration – it does not depend on the
  current parameter estimates).

- cache_dir:

  Directory populated by
  [`trigger_spatial_marks()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_marks.md)
  (built once, reused every iteration – depends only on
  `spatial_trigger`'s grid, which is fixed across iterations).

- cell_size:

  Passed to
  [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md).

- nprocs:

  Passed to
  [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md):
  a whole number of parallel workers (`1`, the default, runs serially),
  `"auto"`, or a cluster from
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).
  A number or `"auto"` starts new workers at every call; to reuse them
  across steps, pass a cluster, as
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  does.

- exclude_self:

  Passed to
  [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md);
  default `TRUE`.

- bgint, trig:

  Optional. The results of
  [`background_intensity()`](https://henriquescadbylab.github.io/spHawkes/reference/background_intensity.md)
  and
  [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md)
  for exactly the `background`, `trend`, `weekly`, `monthly` and
  triggering kernels passed in, as returned in `$bgint` and `$trig` by
  the previous call to `em_step()`. Supplying them saves evaluating the
  same quantities again at the start of this step; the result is the
  same. `NULL` (the default) evaluates them here. Values computed for
  any other inputs give wrong results.

- marks:

  Optional: the spatial marks in `cache_dir`, already read into a list
  with one logical matrix per event, as
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  does once per fit. Passed to
  [`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md);
  `NULL` (the default) reads them from `cache_dir` at this step.

## Value

A list of class `"semihawkes_em_step"` with the updated `mu`, `A`,
`llik`, `lambda_at_all`, `lambda_at_events`, `bgprobs`, the refreshed
`background`, `trend`, `weekly`, `monthly`, `temporal_trigger`,
`spatial_trigger` objects, and `bgint` and `trig`, the background
intensity and triggering contributions evaluated with those refreshed
components – everything needed to feed straight into the next call to
`em_step()`.

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
# One more EM step from where the fit stopped
pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
step <- em_step(fit$events, fit$background, fit$trend, fit$weekly,
                fit$monthly, fit$temporal_trigger, fit$spatial_trigger,
                mu = fit$mu, A = fit$A, pairs = pairs,
                cache_dir = cache_dir, cell_size = 1000)
step
#> <semihawkes_em_step>
#>   mu = 2.55376e-09, A = 0.559732, llik = -1873.06, lambda_at_all = 98
# }
```

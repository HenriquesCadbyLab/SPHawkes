# EM responsibility weights for triggering pairs

Computes `A * g_T(t_i - t_j) * g_S(x_i - x_j, y_i - y_j) / lambda(i)`
for every pair in `pairs` – the numerator of the probability that event
`j` triggered event `i`, i.e. the original script's `excite.wghs`. This
is shared between the temporal and spatial re-estimation steps
([`trigger_reestimate_temporal()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_temporal.md)
/
[`trigger_reestimate_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_spatial.md)),
so it is computed once here rather than duplicated inside each.

## Usage

``` r
trigger_pair_weights(
  events,
  pairs,
  A,
  lambda_at_events,
  temporal_trigger,
  spatial_trigger
)
```

## Arguments

- events:

  A `semihawkes_events` object.

- pairs:

  A pair index from
  [`trigger_pairs()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pairs.md).

- A:

  The current branching-ratio parameter estimate (scalar).

- lambda_at_events:

  The current total-intensity estimate at every event (numeric vector,
  `length(lambda_at_events) == nrow(events$events)` – e.g.
  `background$bg_at_events + A * trigger_evaluate(...)$at_events`).

- temporal_trigger, spatial_trigger:

  The *current* triggering-kernel estimates (from the previous EM
  iteration) used to evaluate the responsibility weights.

## Value

A numeric vector, one weight per row of `pairs`.

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
pairs <- trigger_pairs(fit$events, window_days = 90, window_dist = 2000)
w <- trigger_pair_weights(fit$events, pairs, fit$A, fit$lambda_at_events,
                          fit$temporal_trigger, fit$spatial_trigger)
summary(w)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000000 0.001790 0.006304 0.055504 0.032509 0.979316 
# }
```

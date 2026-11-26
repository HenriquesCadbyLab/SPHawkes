# Re-estimate the temporal triggering kernel from EM responsibility weights

Generalizes the original script's Loop-1 temporal triggering
re-estimation: bins triggering-pair offsets into a weighted histogram
(edge-corrected by
[`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md)'s
output), smooths with
[`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md),
and re-normalizes – returning a `"semihawkes_trigger_temporal"` object,
a drop-in replacement for
[`trigger_temporal_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_temporal_init.md)'s
output (same `base` grid, refreshed `basevalue`/`fun`).

## Usage

``` r
trigger_reestimate_temporal(
  events,
  pairs,
  weights,
  edge_correction,
  temporal_trigger,
  bandwidth = 5
)
```

## Arguments

- events:

  A `semihawkes_events` object.

- pairs:

  A pair index from
  [`trigger_pairs()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pairs.md).

- weights:

  Responsibility weights from
  [`trigger_pair_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pair_weights.md)
  (same length/order as `pairs`).

- edge_correction:

  A `semihawkes_trigger_edge_correction` object from
  [`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md).

- temporal_trigger:

  The current `semihawkes_trigger_temporal` object (only its `base` grid
  is reused; `basevalue`/`fun` are replaced).

- bandwidth:

  Gaussian smoothing bandwidth passed to
  [`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md);
  default `5` (matches the original script).

## Value

A `"semihawkes_trigger_temporal"` object.

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
ec <- trigger_edge_correction(fit$events, fit$temporal_trigger,
                              fit$spatial_trigger, cache_dir)
g_new <- trigger_reestimate_temporal(fit$events, pairs, w, ec,
                                     fit$temporal_trigger)
plot(g_new$base, g_new$basevalue, type = "l", xlim = c(0, 30),
     xlab = "Days since event", ylab = "Density")

# }
```

# Re-estimate the spatial triggering kernel from EM responsibility weights

Spatial triggering re-estimation: bins triggering-pair offsets into a
weighted 2D histogram (edge-corrected by
[`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md)'s
output), smooths with
[`ker.smooth.2D.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.2D.fft.md),
masks to a maximum radius, and re-normalizes – returning a
`"semihawkes_trigger_spatial"` object, a drop-in replacement for
[`trigger_spatial_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_init.md)'s
output.

## Usage

``` r
trigger_reestimate_spatial(
  events,
  pairs,
  weights,
  edge_correction,
  spatial_trigger,
  x_bandwidth = 250,
  y_bandwidth = 250,
  mark_radius = NULL
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
  [`trigger_pair_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pair_weights.md).

- edge_correction:

  A `semihawkes_trigger_edge_correction` object.

- spatial_trigger:

  The current `semihawkes_trigger_spatial` object (its `x`/`y` grid is
  reused; `basevalue`/`fun` are replaced).

- x_bandwidth, y_bandwidth:

  Gaussian smoothing bandwidths passed to
  [`ker.smooth.2D.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.2D.fft.md),
  in the linear units of the working CRS (metres for a UTM projection);
  default `250` for both, i.e. 250 m (matches the original).

- mark_radius:

  Radius beyond which the re-estimated kernel is forced to `0`, in the
  linear units of the working CRS (metres for a UTM projection). The
  default, `NULL`, uses 1.05 times the reach of the kernel's grid (the
  smaller of its largest `x` and `y` offsets), which cuts the kernel to
  a disc just wider than the grid's half-width: 1050 for the default
  +/-1 km grid.

## Value

A `"semihawkes_trigger_spatial"` object.

## Offsets no pair could produce

At grid cells where no event could have produced that spatial offset
while staying inside the region (`repetance == 0`), the weighted density
is zero too, and the ratio of the two is set to `0` rather than `NaN`.
Left as `0 / 0`, a single such cell would spoil the whole kernel through
the FFT smoother.

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
h_new <- trigger_reestimate_spatial(fit$events, pairs, w, ec,
                                    fit$spatial_trigger)
image(h_new$x, h_new$y, h_new$basevalue, asp = 1,
      xlab = "x offset (m)", ylab = "y offset (m)")

# }
```

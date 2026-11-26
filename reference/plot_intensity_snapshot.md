# Plot one gridded intensity snapshot

Plots a previously evaluated background, triggering, or total intensity
snapshot returned by background_trigger_grid().

## Usage

``` r
plot_intensity_snapshot(
  grid,
  time,
  component = c("total", "background", "triggering"),
  cell_size = NULL
)
```

## Arguments

- grid:

  A list returned by background_trigger_grid().

- time:

  A single evaluated model time, in days. This selects a value in
  grid\$days, not a snapshot position.

- component:

  One of "total", "background", or "triggering".

- cell_size:

  Optional positive tile width and height. If NULL, each dimension is
  inferred from the corresponding coordinate spacing. Supply this
  argument when a coordinate axis has only one location.

## Value

A ggplot object.

## Details

Matrix rows correspond to x coordinates and matrix columns correspond to
y coordinates. Matrices are flattened in the same order as expand.grid(x
= grid\$x, y = grid\$y); no transpose is required.

Only grid locations marked as inside the study region are displayed.
Tiles are centred on evaluation coordinates and are not clipped to the
study-region boundary. Values are pointwise intensities, not spatial
integrals over cells.

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
                                days = c(100, 150))
plot_intensity_snapshot(grid, time = 150, component = "total")

# }
```

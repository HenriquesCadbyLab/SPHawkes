# Animate previously evaluated intensity snapshots

Creates a GIF from the snapshots returned by background_trigger_grid().
Each frame displays one evaluated model time. No intermediate intensity
surfaces are generated, and the model is not refitted or reevaluated.

## Usage

``` r
animate_intensity_snapshots(
  grid,
  component = c("total", "background", "triggering"),
  times = NULL,
  filename = tempfile("intensity_", fileext = ".gif"),
  fps = 5,
  width = 900,
  height = 700,
  cell_size = NULL,
  overwrite = FALSE
)
```

## Arguments

- grid:

  A list returned by background_trigger_grid().

- component:

  One of "total", "background", or "triggering".

- times:

  Optional distinct evaluated model times to animate. If NULL, uses all
  times in grid\$days. Selected times are sorted.

- filename:

  Output GIF filename. Defaults to a new file in the R session's
  temporary directory, so nothing is written elsewhere unless you choose
  a location; the path is returned.

- fps:

  Frames per second.

- width, height:

  Output dimensions in pixels.

- cell_size:

  Optional tile size passed to plot_intensity_snapshot().

- overwrite:

  Whether to replace an existing output file.

## Value

The output filename, invisibly.

## Details

All frames share the same colour scale.

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
                                days = seq(100, 130, by = 10))
if (requireNamespace("gifski", quietly = TRUE)) {
  # Written to a temporary file unless `filename` is given
  gif <- animate_intensity_snapshots(grid, component = "triggering",
                                     fps = 2, width = 400, height = 300)
  file.exists(gif)
}
#> [1] TRUE
# }
```

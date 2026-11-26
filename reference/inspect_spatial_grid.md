# Inspect a candidate spatial grid

Construct and optionally plot a regular grid over a spatial window to
inspect a candidate cell size before fitting a model.

## Usage

``` r
inspect_spatial_grid(
  region,
  cell_size,
  points = NULL,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- region:

  An `sf` or `sfc` polygonal spatial window in a Cartesian coordinate
  system with units in metres.

- cell_size:

  A positive number giving the grid cell side length in metres.

- points:

  Optional `sf` or `sfc` point geometries to overlay. Must use the same
  coordinate reference system as `region`.

- plot:

  Logical. Whether to plot the window, grid and points. Default `TRUE`.

- verbose:

  Logical. Whether to print the window dimensions and number of retained
  cells. Default `TRUE`.

## Value

Invisibly returns a list containing:

- grid:

  An `sf` object containing the retained grid cells.

- cell_size:

  The requested cell side length in metres.

- width_km, height_km:

  Width and height of the window's bounding box in kilometres.

- n_cells:

  Number of retained grid cells.

## Details

This function is a descriptive aid, not an automatic cell-size selector.

## Examples

``` r
data("milwaukee")
data("milwaukee_county")

inspection <- inspect_spatial_grid(
  region = milwaukee_county,
  cell_size = 500,
  points = milwaukee
)
#> Window bounding box: 19.84 x 39.07 km
#> Cell size: 500 m
#> Intersecting cells: 2631


inspection$n_cells
#> [1] 2631
```

# Convert a semihawkes background grid to an sf polygon grid

Turns the rectangular
`background$x`/`background$y`/`background$basevalue` grid produced by
[`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md)
into an `sf` object with one square POLYGON per grid cell, so it can be
cropped against `region` with
[`sf::st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
and plotted with
[`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html).

## Usage

``` r
background_to_grid_sf(background, region)
```

## Arguments

- background:

  The `background` component of a `semihawkes_fit` object (has `$x`,
  `$y`, `$basevalue`, `$marks`).

- region:

  An `sf`/`sfc` POLYGON/MULTIPOLYGON, used only to supply the CRS.

## Value

An `sf` object with columns `x`, `y`, `basevalue`, and one square
POLYGON geometry per grid cell.

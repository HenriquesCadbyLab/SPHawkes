# Crop a semihawkes spatial background to its study-area sf polygon

Restricts the fitted spatial-background grid to the `region` the model
was fitted on, instead of the full rectangular estimation grid, by
intersecting each grid cell (as a small square polygon) with `region`.

## Usage

``` r
crop_background_spatial(background, region, crop = TRUE)
```

## Arguments

- background:

  The `background` component of a `semihawkes_fit` object.

- region:

  An `sf`/`sfc` POLYGON/MULTIPOLYGON study-area boundary, typically
  `fit$events$region`.

- crop:

  Currently unused; kept for interface compatibility. Cells outside
  `region` are always dropped (they have no intersection).

## Value

An `sf` object with columns `x`, `y`, `basevalue`, and one (possibly
clipped) polygon geometry per grid cell that overlaps `region`, ready
for `ggplot2::geom_sf(aes(fill = basevalue))`.

## Details

[`sf::st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
between polygons can occasionally return mixed geometry types (e.g. a
`POLYGON` clipped down to a `LINESTRING` or `POINT` when it only grazes
the boundary, or a `GEOMETRYCOLLECTION`).
[`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html)/`fortify()`
cannot handle a geometry column with mixed/non-standard types, so we
explicitly keep only (MULTI)POLYGON results and cast them to
POLYGON/MULTIPOLYGON.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
bg <- background_spatial(ev, cell_size = 1000)
cells <- crop_background_spatial(bg, ev$region)
head(cells)
#> Simple feature collection with 6 features and 3 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 5500 ymax: 500
#> Projected CRS: WGS 84 / UTM zone 17S
#>      x y basevalue                       geometry
#> 1    0 0 0.1085815 MULTIPOLYGON (((500 500, 50...
#> 2 1000 0 0.1176292 MULTIPOLYGON (((500 500, 15...
#> 3 2000 0 0.1557842 MULTIPOLYGON (((1500 500, 2...
#> 4 3000 0 0.2630437 MULTIPOLYGON (((2500 500, 3...
#> 5 4000 0 0.3389455 MULTIPOLYGON (((3500 500, 4...
#> 6 5000 0 0.3366383 MULTIPOLYGON (((4500 500, 5...
```

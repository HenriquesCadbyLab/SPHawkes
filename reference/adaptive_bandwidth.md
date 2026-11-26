# Adaptive (nearest-neighbour) spatial bandwidth for each event

For every point, the bandwidth is the distance to its `k`-th nearest
neighbour among neighbours at least `min_dist` away (this skips
near-duplicate/co-located points, e.g. events geocoded to the same
address), floored at `min_bandwidth`.

## Usage

``` r
adaptive_bandwidth(
  x,
  y,
  k = 5,
  min_dist = 150,
  min_bandwidth = 250,
  k_search = 50
)
```

## Arguments

- x, y:

  Numeric vectors of projected event coordinates (same length).

- k:

  Which nearest neighbour (beyond `min_dist`) to use; default `5`.

- min_dist:

  Minimum distance for a neighbour to be considered (neighbours closer
  than this are ignored), in the linear units of the working CRS (metres
  for a UTM projection); default `150`, i.e. 150 m.

- min_bandwidth:

  Floor applied to every bandwidth, in the linear units of the working
  CRS (metres for a UTM projection); default `250`, i.e. 250 m.

- k_search:

  How many nearest neighbours to search over before filtering by
  `min_dist`; increase this if a warning about insufficient neighbours
  is raised. Default `50`.

## Value

A numeric vector of bandwidths, one per input point.

## Examples

``` r
xy <- sf::st_coordinates(sim_events)
bw <- adaptive_bandwidth(xy[, 1], xy[, 2], k = 5, min_dist = 150,
                         min_bandwidth = 300)
summary(bw)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   300.0   300.0   352.3   469.8   507.7  2057.0 
```

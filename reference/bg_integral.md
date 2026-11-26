# Monte Carlo integral of each event's isotropic Gaussian kernel over the study region

For every event, its background kernel (an isotropic Gaussian centred on
the event with standard deviation `bandwidth[i]`) generally has some of
its mass falling outside the study region boundary. This computes the
fraction of that mass that falls *inside* `region`, by Monte Carlo
integration over `n_mc` random points drawn uniformly in the region. It
is used to edge-correct the background-density kernel estimate.

## Usage

``` r
bg_integral(
  coorx,
  coory,
  bandwidth,
  region,
  n_mc = 3000,
  cap = 0.99,
  seed = NULL,
  restore_rng = TRUE
)
```

## Arguments

- coorx, coory:

  Numeric vectors of projected event coordinates.

- bandwidth:

  Numeric vector (same length), the per-event standard deviation from
  e.g.
  [`adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/adaptive_bandwidth.md).

- region:

  An `sf`/`sfc` polygon (already in the same projected CRS as the
  coordinates).

- n_mc:

  Number of Monte Carlo points drawn in `region`; default `3000`,
  matching the original script.

- cap:

  Values above `cap` are truncated to `1` (default `0.99`), since a
  fully-interior kernel should integrate to ~1 and larger Monte Carlo
  noise is expected near that ceiling.

- seed:

  Optional random seed, making the Monte Carlo draw reproducible.

- restore_rng:

  Only used when `seed` is supplied. If `TRUE` (the default), the global
  random-number state is put back as it was before the call, so the seed
  affects only this call. If `FALSE`, the call leaves the generator in
  the state `set.seed(seed)` and its own draws produce, as earlier
  versions of the package did – useful for reproducing results from
  scripts that relied on that.

## Value

A numeric vector (same length as `coorx`) with each event's in-region
kernel mass, in `(0, 1]`.

## Examples

``` r
xy <- sf::st_coordinates(sim_events)[1:5, ]
# Share of a 300 m Gaussian kernel around each event that lies in the region
bg_integral(xy[, 1], xy[, 2], bandwidth = rep(300, 5), region = sim_region,
            n_mc = 500, seed = 1)
#> [1] 0.6363517 1.0000000 1.0000000 1.0000000 0.5050622
```

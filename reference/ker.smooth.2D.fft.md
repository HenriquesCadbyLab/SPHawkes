# 2D Gaussian kernel smoothing by FFT, with reflected boundaries

Smooths a matrix `z` sampled on an equally-spaced `x`/`y` grid with a
separable Gaussian kernel of standard deviations
`x.bandwidth`/`y.bandwidth`, using an FFT.

## Usage

``` r
ker.smooth.2D.fft(x, y, z, x.bandwidth, y.bandwidth)
```

## Arguments

- x, y:

  Numeric vectors of equally-spaced grid points for each dimension.

- z:

  Numeric matrix of values on the `x` by `y` grid
  (`nrow(z) == length(x)`, `ncol(z) == length(y)`).

- x.bandwidth, y.bandwidth:

  Standard deviations of the Gaussian smoothing kernel along `x` and
  `y`.

## Value

A numeric matrix of dimension `(length(x)+1) x (length(y)+1)`: the
smoothed values at the points halfway between successive `x` (and `y`),
and half a grid step beyond the first and last. When `x` and `y` are the
midpoints of histogram bins, as from
[`hist.weighted.2D()`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.2D.md),
these are the bin edges.

## Details

Note that despite the shared "fft" in the name, this extends `z` by
*reflection* in both dimensions, whereas the 1D
[`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md)
extends periodically — appropriately so, since this one smooths the
spatial triggering kernel (not periodic) and that one smooths the weekly
and monthly background components (periodic). See
[`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md).

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
# Smooth a single spike at the origin; the result, one more row and column
# than the input, is centred on the spike
g <- seq(-1000, 1000, by = 100)
z <- matrix(0, length(g), length(g))
z[11, 11] <- 1
smooth <- ker.smooth.2D.fft(g, g, z, x.bandwidth = 200, y.bandwidth = 200)
dim(smooth)
#> [1] 22 22
# The largest values sit half a step either side of the spike
edges <- c(g - 50, max(g) + 50)
peak <- which(smooth > max(smooth) * (1 - 1e-9), arr.ind = TRUE)
edges[unique(peak[, 1])]
#> [1] -50  50
```

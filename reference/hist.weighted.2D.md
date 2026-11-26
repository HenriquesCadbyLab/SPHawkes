# Weighted 2D histogram

2D analogue of
[`hist.weighted()`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md):
bins `(x, y)` pairs into a grid defined by `x.breaks`/`y.breaks` and
sums `weight` within each cell.

## Usage

``` r
# S3 method for class 'weighted.2D'
hist(x, y, weight, x.breaks, y.breaks)
```

## Arguments

- x, y:

  Numeric vectors of coordinates to bin (same length).

- weight:

  Numeric vector of weights, same length as `x`/`y`.

- x.breaks, y.breaks:

  Numeric vectors of bin edges for each dimension; must span the range
  of `x`/`y` respectively (an error otherwise).

## Value

A list with `x.breaks`, `y.breaks`, `x.mids`, `y.mids`, `counts` (matrix
of summed weight per cell), and `density` (`counts` normalised to
integrate to 1 over the grid).

## A note on the name

As with
[`hist.weighted()`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md),
this is an ordinary function rather than an S3 method for
[`graphics::hist()`](https://rdrr.io/r/graphics/hist.html), and needs an
explicit `@rawNamespace` export directive for that reason.

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
h <- hist.weighted.2D(
  x = c(0.2, 0.7, 1.5), y = c(0.1, 0.9, 1.2), weight = c(1, 2, 1),
  x.breaks = 0:2, y.breaks = 0:2
)
h$density
#>      [,1] [,2]
#> [1,] 0.75 0.00
#> [2,] 0.00 0.25
```

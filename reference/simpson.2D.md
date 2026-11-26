# Simpson's rule double integral of an equally-spaced function

2D analogue of
[`simpson()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.md):
approximates \\\int\int f(x,y)\\dx\\dy\\ from a matrix of function
values sampled on an equally-spaced grid in both dimensions.

## Usage

``` r
simpson.2D(value, dx = 1, dy = 1)
```

## Arguments

- value:

  Numeric matrix of function values, rows indexing the first dimension
  (step `dx`) and columns the second (step `dy`). As for
  [`simpson()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.md),
  both the number of rows and the number of columns must be odd and at
  least 3.

- dx, dy:

  Grid spacing along rows/columns; default `1` each.

## Value

A single numeric value, the approximate double integral.

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
g <- seq(-4, 4, by = 0.1)
z <- outer(stats::dnorm(g), stats::dnorm(g))
simpson.2D(z, dx = 0.1, dy = 0.1)
#> [1] 0.9998733
```

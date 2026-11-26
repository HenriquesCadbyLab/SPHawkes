# Simpson's rule integral of an equally-spaced function

Approximates \\\int f(x)\\dx\\ from function values `value` sampled on
an equally-spaced grid with step `dx`, using (composite) Simpson's rule.

## Usage

``` r
simpson(value, dx = 1)
```

## Arguments

- value:

  Numeric vector of function values on an equally-spaced grid. Simpson's
  rule needs an odd number of points (an even number of intervals), at
  least 3; any other length is an error.

- dx:

  Grid spacing; default `1`.

## Value

A single numeric value, the approximate integral.

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
# The standard normal density integrates to (almost exactly) 1
x <- seq(-4, 4, by = 0.1)
simpson(stats::dnorm(x), dx = 0.1)
#> [1] 0.9999366
```

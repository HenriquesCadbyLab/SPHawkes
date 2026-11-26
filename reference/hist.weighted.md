# Weighted 1D histogram

Bins `x` into `breaks` and sums `weight` (instead of counting) within
each bin, returning both counts and a weighted density (counts divided
by total weight and by bin width) plus bin midpoints — the building
block for the non-parametric weekly/monthly/trend estimators.

## Usage

``` r
# S3 method for class 'weighted'
hist(x, weight, breaks)

# S3 method for class 'weighted.1D'
hist(x, weight, breaks)
```

## Arguments

- x:

  Numeric vector of values to bin.

- weight:

  Numeric vector of weights, same length as `x` (e.g. the EM
  background/trigger probabilities).

- breaks:

  Numeric vector of bin edges; must span the range of `x` (an error
  otherwise).

## Value

A list with `breaks`, `counts` (summed weight per bin), `density`
(`counts` normalised to integrate to 1 over `breaks`), and `mids` (bin
midpoints).

## A note on the name

Despite the dots, this is an ordinary function, **not** an S3 method for
the [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html) generic —
the name is inherited from the original D'Angelo et al. (2024) script.
It has to be exported with an explicit `@rawNamespace` directive because
roxygen2 otherwise reads `hist.weighted` as `hist` dispatching on class
`"weighted"` and registers it as a method, which both hides the function
from users and makes [`hist()`](https://rdrr.io/r/graphics/hist.html)
misbehave on any object that happens to carry that class.

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
x <- c(0.2, 0.4, 1.3, 2.7, 2.9)
w <- c(1, 2, 1, 1, 3)
h <- hist.weighted(x, w, breaks = 0:3)
h$counts
#> 1 2 3 
#> 3 1 4 
h$density
#>     1     2     3 
#> 0.375 0.125 0.500 
```

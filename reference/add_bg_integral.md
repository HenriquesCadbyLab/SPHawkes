# Add the background-kernel edge-correction integral to a `semihawkes_events` object

Convenience wrapper around
[`bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/bg_integral.md).
Requires `x$events$bandwidth` to already exist (see
[`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md)),
and additionally drops any event that falls outside `x$region` after
re-checking containment.

## Usage

``` r
add_bg_integral(x, ...)
```

## Arguments

- x:

  A `semihawkes_events` object with a `bandwidth` column already added.

- ...:

  Passed on to
  [`bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/bg_integral.md)
  (`n_mc`, `cap`, `seed`, `restore_rng`).

## Value

`x`, with a `bg.integral` column added to `x$events`.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
summary(ev$events$bg.integral)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2892  0.6479  0.8855  0.8137  1.0000  1.0000 
```

# Add an adaptive bandwidth column to a `semihawkes_events` object

Convenience wrapper around
[`adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/adaptive_bandwidth.md)
that operates directly on the `events` data.frame produced by
[`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md).

## Usage

``` r
add_adaptive_bandwidth(x, ...)
```

## Arguments

- x:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)).

- ...:

  Passed on to
  [`adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/adaptive_bandwidth.md)
  (`k`, `min_dist`, `min_bandwidth`, `k_search`).

## Value

`x`, with a `bandwidth` column added to `x$events`.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
summary(ev$events$bandwidth)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   300.0   300.0   352.3   469.8   507.7  2057.0 
```

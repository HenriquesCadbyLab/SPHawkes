# Build the (i, j) triggering-pair index for re-estimation

Enumerates ordered event pairs `(i, j)` where `j` could plausibly have
triggered `i`: `j` happened strictly before `i`, within `window_days`,
and within `window_dist` in both spatial dimensions.

## Usage

``` r
trigger_pairs(events, window_days = 90, window_dist = 1000)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)).

- window_days:

  Maximum number of days after `j` that `i` can occur; default `90`
  (matches the original script's actual filter, `<= 90.0`).

- window_dist:

  Maximum absolute distance (in both x and y) between `i` and `j`, in
  the linear units of the working CRS (metres for a UTM projection);
  default `1000`, i.e. 1 km, the same as the default reach of the
  spatial triggering kernel (see *Choosing the triggering reach* in
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)).

## Value

A data frame with integer columns `i`, `j` (row indices into
`events$events`), one row per qualifying pair.

## Scalability

Once the events are sorted by day, the events that `j` could have
triggered (those after it and within `window_days`) form one run, found
with [`findInterval()`](https://rdrr.io/r/base/findInterval.html); only
those are checked against `window_dist`. Time and memory therefore grow
with the number of candidate pairs rather than with `n^2`, as they did
when every pair was compared. Pairs are returned in the same order as
before: by `j`, then by `i`.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
pairs <- trigger_pairs(ev, window_days = 30, window_dist = 1000)
nrow(pairs)
#> [1] 667
head(pairs)
#>   i j
#> 1 6 1
#> 2 8 1
#> 3 3 2
#> 4 4 2
#> 5 5 2
#> 6 7 2
```

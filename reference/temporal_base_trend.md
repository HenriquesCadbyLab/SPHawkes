# Long-term trend temporal base value

Estimates a non-periodic, long-term background multiplier over the whole
study horizon. With `resp_weight = NULL` (default) this returns a flat
trend (`basevalue = 1` everywhere) - the Loop-0 initial guess in the
original script, where the weighted-smoothing loop was disabled and the
trend was left flat until the first E-step. With `resp_weight` supplied
(one EM responsibility weight per event, e.g.
`trend.fun(days) * mub.events / lambda.at.events`), it reproduces the
Loop-k weighted Gaussian smoother with boundary (edge-effect) correction
via [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html).

## Usage

``` r
temporal_base_trend(
  events,
  resp_weight = NULL,
  bandwidth = 150,
  chunk_size = 500
)
```

## Arguments

- events:

  A `"semihawkes_events"` object from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md).

- resp_weight:

  Optional numeric vector, one entry per event in `events$events`, e.g.
  an EM E-step responsibility weight. If `NULL` (default), only the
  day-count correction weight is used (the Loop-0 initial-guess
  behaviour).

- bandwidth:

  Gaussian smoothing bandwidth (in days); default `150`, matching the
  original script.

- chunk_size:

  Number of distinct event days processed per chunk when `resp_weight`
  is supplied, to bound peak memory use for long study periods; default
  `500`.

## Value

A list of class `"semihawkes_temporal_base"` with elements `breaks`
(`events$time.marks`), `basevalue` (normalized to mean 1), `period`
(`NA`, i.e. non-periodic), and `fun` (a vectorized function of `days`,
`0` outside `[0, TT]`).

## Details

Events fall on whole days, so the weights of events on the same day are
added up first, and the smoothed curve is evaluated on a grid of one
day, or `bandwidth / 100` days if that is finer, then interpolated
linearly to `events$time.marks`. That changes the curve by at most about
`(step / bandwidth)^2 / 8` of its value: 6 parts in a million for the
default bandwidth. When the grid would be as fine as the time marks, the
curve is evaluated at the time marks directly.

## References

Zhuang & Mateu (2019); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
# Without weights: the flat initial guess used before the first EM step
flat <- temporal_base_trend(ev)
range(flat$basevalue)
#> [1] 1 1
# With weights (here every event weighted equally): a smoothed trend
trend <- temporal_base_trend(ev, resp_weight = rep(1, nrow(ev$events)))
plot(trend$breaks, trend$basevalue, type = "l",
     xlab = "Days since start", ylab = "Trend multiplier")
```

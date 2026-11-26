# Weekly non-parametric temporal base value

Estimates a periodic (7-day) background multiplier from event
days-of-week, following the weekly block of the original Zhuang & Mateu
(2019) style script: bin events by day-of-week with
[`hist.weighted()`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md),
smooth with
[`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md),
and normalize to a mean of 1. Generalizes both the Loop-0 initial guess
(`resp_weight = NULL`) and the Loop-k E-step re-estimation
(`resp_weight` = the EM responsibility weight for the weekly component
at each event, e.g. `weekly.fun(days) * mub.events / lambda.at.events`
in the original script) into one function.

## Usage

``` r
temporal_base_weekly(
  events,
  resp_weight = NULL,
  breaks = seq(0, 7, 0.05),
  bandwidth = 1
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

- breaks:

  Numeric grid of day-of-week bin edges; default `seq(0, 7, 0.05)` (`0`
  = Sunday), matching the original script.

- bandwidth:

  Gaussian smoothing bandwidth passed to
  [`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md);
  default `1`.

## Value

A list of class `"semihawkes_temporal_base"` with elements `breaks`,
`basevalue` (normalized to mean 1), `period` (`7`), and `fun` (a
vectorized function of `days` returning the periodic multiplier).

## References

Zhuang & Mateu (2019); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
weekly <- temporal_base_weekly(ev)
weekly
#> <semihawkes_temporal_base>
#>   Grid points : 141 (range 0 to 7)
#>   Period      : 7
#>   Base value  : mean 1, range [0.839, 1.2]
# Multiplier for each day of the week (0 = Sunday)
round(weekly$fun(0:6), 2)
#> [1] 1.19 1.04 0.88 0.85 0.91 0.99 1.13
```

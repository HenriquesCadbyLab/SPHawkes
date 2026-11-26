# Monthly non-parametric temporal base value

Estimates a periodic (12-month) background multiplier from each event's
calendar month (`events$events$month`, `0` = January), the same way as
[`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md)
but with a 12-month period. See that function's documentation for the
`resp_weight` argument and the Loop-0 vs. Loop-k generalization.

## Usage

``` r
temporal_base_monthly(
  events,
  resp_weight = NULL,
  breaks = seq(0, 12, 0.05),
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

  Numeric grid of month-of-year bin edges; default `seq(0, 12, 0.05)`.

- bandwidth:

  Gaussian smoothing bandwidth passed to
  [`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md);
  default `1`.

## Value

A list of class `"semihawkes_temporal_base"`; see
[`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md).

## References

Zhuang & Mateu (2019); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
monthly <- temporal_base_monthly(ev)
# Multiplier for each calendar month (0 = January)
round(monthly$fun(0:11), 2)
#>  [1] 0.88 0.86 1.03 1.14 1.04 0.96 0.97 1.05 1.08 1.04 1.00 0.96
```

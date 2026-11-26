# Initial (parametric) temporal triggering-function guess

Builds the Loop-0 power-law temporal triggering kernel \\g(t) =
(\text{offset} + t/\text{scale})^{\text{exponent}}\\, normalized to
integrate to `1` over `[0, cutoff]` via
[`simpson()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.md).

## Usage

``` r
trigger_temporal_init(
  cutoff = 90,
  step = 0.05,
  offset = 0.1,
  scale = 10,
  exponent = -1.03
)
```

## Arguments

- cutoff:

  Maximum time lag considered, in days; default `90`.

- step:

  Grid step for `base`, in days; default `0.05`. `cutoff / step` must be
  an even whole number (1800 for the defaults), so that the grid ends at
  `cutoff` and has the odd number of points
  [`simpson()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.md)
  needs.

- offset, scale, exponent:

  Power-law shape parameters: \\g(t) = (\text{offset} +
  t/\text{scale})^{\text{exponent}}\\. Defaults (`0.1`, `10`, `-1.03`)
  match the original script.

## Value

A list of class `"semihawkes_trigger_temporal"` with `base` (the
time-lag grid), `basevalue` (normalized kernel values), and `fun` (a
vectorized evaluator, `0` outside `[0, cutoff]`).

## Examples

``` r
g <- trigger_temporal_init()
g
#> <semihawkes_trigger_temporal>
#>   1801 points over [0.00, 90.00] days
plot(g$base, g$basevalue, type = "l", log = "y",
     xlab = "Days since event", ylab = "Density")
```

# 1D Gaussian kernel smoothing by FFT, with periodic boundaries

Smooths values `z` sampled at equally-spaced points `x` with a Gaussian
kernel of standard deviation `bandwidth`, using an FFT.

## Usage

``` r
ker.smooth.fft(x, z, bandwidth)

ker.smooth.1D.fft(x, z, bandwidth)
```

## Arguments

- x:

  Numeric vector of equally-spaced grid points.

- z:

  Numeric vector of values at `x` (same length).

- bandwidth:

  Standard deviation of the Gaussian smoothing kernel.

## Value

A numeric vector of length `length(z) + 1`: the smoothed values at the
points halfway between successive `x`, and half a grid step before the
first and after the last. When `x` are the midpoints of histogram bins,
as from
[`hist.weighted()`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md),
these are the bin edges.

## Not interchangeable with [`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md)

This is **not** simply a faster version of
[`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md):
the two treat the boundaries differently, and generally return different
values.

- [`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md)
  extends the series by *reflection* (`c(y, rev(y))`), which is the
  right assumption for a curve that is not periodic — such as the
  temporal triggering kernel, whose lag axis simply stops at the cutoff.

- `ker.smooth.fft()` extends the series *periodically* (`c(z, z)`),
  which is the right assumption for a genuinely cyclical curve — such as
  the day-of-week and month-of-year background components, where the
  value at Saturday midnight really does wrap round to Sunday.

The package uses each accordingly:
[`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md)
and
[`temporal_base_monthly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_monthly.md)
call this function, while
[`trigger_reestimate_temporal()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_temporal.md)
calls
[`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md),
and
[`trigger_reestimate_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_spatial.md)
calls
[`ker.smooth.2D.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.2D.fft.md)
(which, despite also being FFT-based, reflects rather than wraps).
Substituting one for another "for speed" will silently change the
estimates near the boundaries.

## References

D'Angelo et al. (2024); see
[spHawkes-package](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md).

## Examples

``` r
# Periodic smoothing over a week; the result has one value per break
breaks <- seq(0, 7, by = 0.05)
h <- hist.weighted(c(1, 1, 2, 5, 5, 5), rep(1, 6), breaks = breaks)
smooth <- ker.smooth.fft(h$mids, h$density, bandwidth = 0.5)
plot(breaks, smooth, type = "l", xlab = "Day of week", ylab = "Smoothed density")
```

# Check consistency of background probabilities and excitation weights

For each event `i`, verifies that
`bgprobs[i] + sum(pair_weights[pairs$i == i])` equals `1` (up to a
numerical tolerance) – the condition
[`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)
relies on for `u ~ Uniform(0, 1)` to be partitioned correctly between
the background outcome and every candidate-parent outcome. Call this
once before drawing samples with
[`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)/[`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md).

## Usage

``` r
check_branching_weights(bgprobs, pairs, pair_weights, tol = 1e-06)
```

## Arguments

- bgprobs, pairs, pair_weights:

  Same as in
  [`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md).

- tol:

  Numerical tolerance for the sum-to-1 check; default `1e-6`.

## Value

A list of class `"semihawkes_branching_check"` with elements:

- totals:

  Numeric vector, one entry per event:
  `bgprobs[i] + sum(pair_weights[pairs$i == i])`.

- max_abs_error:

  `max(abs(totals - 1))`.

- n_no_parents:

  Number of events with `bgprobs[i] < 1` but no row in `pairs`
  (candidate parents) – these silently fall back to background in
  [`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md),
  which is only correct if `bgprobs[i]` is actually `1` for those
  events.

- ok:

  `TRUE` if `max_abs_error <= tol` and `n_no_parents == 0`.

## Examples

``` r
bgprobs <- c(1, 0.4, 0.5, 0.2)
pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
                    j = as.integer(c(1, 1, 2, 2, 3)))
pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)

check <- check_branching_weights(bgprobs, pairs, pair_weights)
check
#> <semihawkes_branching_check>
#>   max |total - 1| = 0
#>   events with no candidate parents (bgprobs < 1): 0
#>   ok: TRUE
stopifnot(check$ok)
```

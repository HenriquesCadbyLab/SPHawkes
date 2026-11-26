# Draw multiple independent realizations of the branching structure

Repeatedly calls
[`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)
to build a Monte Carlo ensemble of branching-structure realizations, for
quantifying uncertainty in any summary that depends on the
background/triggered classification of events (e.g. the number of
background events per month, or their spatial concentration – see the
package vignette on branching-structure UQ).

## Usage

``` r
sample_branching_multi(
  bgprobs,
  pairs,
  pair_weights,
  n_samples = 100,
  seed = NULL,
  restore_rng = TRUE
)
```

## Arguments

- bgprobs, pairs, pair_weights:

  Passed to
  [`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)
  on every draw (fixed across draws – only the internal
  `u ~ Uniform(0,1)` draws differ between samples).

- n_samples:

  Number of independent realizations to draw; default `100`.

- seed:

  Optional seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing, for
  reproducibility. Default `NULL` (no seed set – uses the current RNG
  state).

- restore_rng:

  Only used when `seed` is supplied. If `TRUE` (the default), the global
  random-number state is put back as it was before the call, so the seed
  affects only this call. If `FALSE`, the call leaves the generator in
  the state `set.seed(seed)` and its own draws produce, as earlier
  versions of the package did – useful for reproducing results from
  scripts that relied on that.

## Value

A list of length `n_samples`, each element a
`"semihawkes_branching_sample"` object from
[`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md).

## See also

[`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)
for a single draw,
[`summarize_branching_samples()`](https://henriquescadbylab.github.io/spHawkes/reference/summarize_branching_samples.md)
to extract per-sample summaries (e.g. number of background events) from
the ensemble.

## Examples

``` r
bgprobs <- c(1, 0.4, 0.5, 0.2)
pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
                    j = as.integer(c(1, 1, 2, 2, 3)))
pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)

samples <- sample_branching_multi(bgprobs, pairs, pair_weights,
                                  n_samples = 200, seed = 1)
length(samples)
#> [1] 200
summarize_branching_samples(samples)
#> <semihawkes_branching_summary>
#>   200 samples
#>   median background events: 2
#>   IQR: [2, 3]
```

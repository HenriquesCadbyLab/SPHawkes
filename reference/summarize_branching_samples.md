# Summarize an ensemble of branching-structure samples

Extracts the number of background events from each realization in a
[`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md)
ensemble, for reporting the central tendency of the background-event
count under sampling uncertainty. Reports the median only (not the mean
or any fixed-threshold count such as `sum(bgprobs >= 0.5)`), since the
empirical distribution of background counts is often asymmetric and the
median is more robust to that asymmetry.

## Usage

``` r
summarize_branching_samples(samples)
```

## Arguments

- samples:

  A list of `"semihawkes_branching_sample"` objects, from
  [`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md).

## Value

A list of class `"semihawkes_branching_summary"` with elements:

- num_bg:

  Integer vector, length `length(samples)`: the number of background
  events in each sample.

- median_bg:

  `median(num_bg)`.

- iqr_bg:

  Named numeric vector, `quantile(num_bg, c(0.25, 0.75))`.

## Examples

``` r
bgprobs <- c(1, 0.4, 0.5, 0.2)
pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
                    j = as.integer(c(1, 1, 2, 2, 3)))
pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)

samples <- sample_branching_multi(bgprobs, pairs, pair_weights,
                                  n_samples = 200, seed = 1)
summarize_branching_samples(samples)
#> <semihawkes_branching_summary>
#>   200 samples
#>   median background events: 2
#>   IQR: [2, 3]
```

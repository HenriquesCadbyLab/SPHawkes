# Sample one realization of the branching structure via UQ

Draws a single stochastic realization of the Hawkes branching structure
by, for each event `i`, sampling whether it is a background event or was
triggered by a specific earlier event `j`, using the fitted background
probabilities and pairwise excitation weights. This is the "stochastic
declustering" step: unlike thresholding `bgprobs` at a fixed cutoff
(e.g. `0.5`), it produces a proper Monte Carlo draw from the branching
structure's posterior distribution given the fitted model, so repeated
calls can be used to quantify uncertainty in any downstream summary that
depends on which events are "background" vs. "triggered" (see
[`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md)).

## Usage

``` r
sample_branching(bgprobs, pairs, pair_weights)
```

## Arguments

- bgprobs:

  Numeric vector of fitted background probabilities, one per event, in
  the same row order as `events$events` (e.g. `fit$bgprobs` from
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)/[`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md)).

- pairs:

  A pair index from
  [`trigger_pairs()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pairs.md),
  with integer columns `i` (child/triggered event) and `j` (candidate
  parent event). Must be the *same* `pairs` object used to fit the fit
  whose `bgprobs` and `pair_weights` are supplied – the window
  parameters (`window_days`/`window_dist`) determine which parents are
  even candidates.

- pair_weights:

  Numeric vector of excitation weights from
  [`trigger_pair_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pair_weights.md),
  same length and row order as `pairs`. For each child `i`,
  `bgprobs[i] + sum(pair_weights[pairs$i == i])` should equal `1` (up to
  floating-point error); see
  [`check_branching_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/check_branching_weights.md).

## Value

A list of class `"semihawkes_branching_sample"` with a single element
`mark`: an integer vector, one entry per event, `0` if that event was
sampled as background, or the row index (into `events$events`) of its
sampled parent otherwise.

## Details

For each event `i`:

- Draw `u ~ Uniform(0, 1)`.

- If `u < bgprobs[i]`, mark `i` as background (`mark[i] = 0`).

- Otherwise, find all candidate parents `j` of `i` (rows of `pairs` with
  `pairs$i == i`), and assign `i` to the parent `j` whose cumulative
  excitation-weight interval contains `u`, i.e. sample `j` with
  probability proportional to `pair_weights` among `i`'s candidate
  parents.

## See also

[`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md)
for repeated draws,
[`check_branching_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/check_branching_weights.md)
to validate `bgprobs`/`pair_weights` consistency before sampling.

## Examples

``` r
# A tiny synthetic case: 4 events, with `bgprobs` and `pair_weights`
# chosen so each event's total probability sums to 1 (see
# `check_branching_weights()`). Event 1 has no candidate parents, so it
# is always background.
bgprobs <- c(1, 0.4, 0.5, 0.2)
pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
                    j = as.integer(c(1, 1, 2, 2, 3)))
pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)

set.seed(1)
draw <- sample_branching(bgprobs, pairs, pair_weights)
draw
#> <semihawkes_branching_sample>
#>   4 events, 2 background (50.0%), 2 triggered
draw$mark # 0 = background; otherwise the sampled parent's row index
#> [1] 0 0 1 3

# See `vignette("branching-structure-uq")` for the full workflow on a
# fitted `semihawkes_fit` object.
```

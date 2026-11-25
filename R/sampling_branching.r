#' Sample one realization of the branching structure via UQ
#'
#' Draws a single stochastic realization of the Hawkes branching structure
#' by, for each event `i`, sampling whether it is a background event or was
#' triggered by a specific earlier event `j`, using the fitted background
#' probabilities and pairwise excitation weights. This is the "stochastic
#' declustering" step: unlike thresholding `bgprobs` at a fixed cutoff
#' (e.g. `0.5`), it produces a proper Monte Carlo draw from the branching
#' structure's posterior distribution given the fitted model, so repeated
#' calls can be used to quantify uncertainty in any downstream summary that
#' depends on which events are "background" vs. "triggered" (see
#' `sample_branching_multi()`).
#'
#' For each event `i`:
#' \itemize{
#'   \item Draw `u ~ Uniform(0, 1)`.
#'   \item If `u < bgprobs[i]`, mark `i` as background (`mark[i] = 0`).
#'   \item Otherwise, find all candidate parents `j` of `i` (rows of
#'     `pairs` with `pairs$i == i`), and assign `i` to the parent `j` whose
#'     cumulative excitation-weight interval contains `u`, i.e. sample `j`
#'     with probability proportional to `pair_weights` among `i`'s
#'     candidate parents.
#' }
#'
#' @param bgprobs Numeric vector of fitted background probabilities, one
#'   per event, in the same row order as `events$events` (e.g.
#'   `fit$bgprobs` from `fit_semihawkes()`/`em_step()`).
#' @param pairs A pair index from `trigger_pairs()`, with integer columns
#'   `i` (child/triggered event) and `j` (candidate parent event). Must be
#'   the *same* `pairs` object used to fit the fit whose `bgprobs` and
#'   `pair_weights` are supplied -- the window parameters
#'   (`window_days`/`window_dist`) determine which parents are even
#'   candidates.
#' @param pair_weights Numeric vector of excitation weights from
#'   `trigger_pair_weights()`, same length and row order as `pairs`. For
#'   each child `i`, `bgprobs[i] + sum(pair_weights[pairs$i == i])` should
#'   equal `1` (up to floating-point error); see
#'   `check_branching_weights()`.
#'
#' @return A list of class `"semihawkes_branching_sample"` with a single
#'   element `mark`: an integer vector, one entry per event, `0` if that
#'   event was sampled as background, or the row index (into
#'   `events$events`) of its sampled parent otherwise.
#'
#' @seealso [sample_branching_multi()] for repeated draws,
#'   [check_branching_weights()] to validate `bgprobs`/`pair_weights`
#'   consistency before sampling.
#'
#' @examples
#' # A tiny synthetic case: 4 events, with `bgprobs` and `pair_weights`
#' # chosen so each event's total probability sums to 1 (see
#' # `check_branching_weights()`). Event 1 has no candidate parents, so it
#' # is always background.
#' bgprobs <- c(1, 0.4, 0.5, 0.2)
#' pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
#'                     j = as.integer(c(1, 1, 2, 2, 3)))
#' pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)
#'
#' set.seed(1)
#' draw <- sample_branching(bgprobs, pairs, pair_weights)
#' draw
#' draw$mark # 0 = background; otherwise the sampled parent's row index
#'
#' # See `vignette("branching-structure-uq")` for the full workflow on a
#' # fitted `semihawkes_fit` object.
#' @export
sample_branching <- function(bgprobs, pairs, pair_weights) {
  if (!is.numeric(bgprobs) || any(!is.finite(bgprobs)) || anyNA(bgprobs)) {
    stop("`bgprobs` must be a numeric vector.", call. = FALSE)
  }
  if (any(bgprobs < 0 | bgprobs > 1)) {
    stop("`bgprobs` must contain probabilities in [0, 1].", call. = FALSE)
  }
  if (!is.data.frame(pairs) ||
      !all(c("i", "j") %in% names(pairs)) ||
      !is.integer(pairs$i) || !is.integer(pairs$j)) {
    stop("`pairs` must be a data frame with integer columns `i` and `j` ",
         "(see `trigger_pairs()`).", call. = FALSE)
  }
  if (!is.numeric(pair_weights) || length(pair_weights) != nrow(pairs) || anyNA(pair_weights)) {
    stop("`pair_weights` must be numeric with length equal to `nrow(pairs)` ",
         "(see `trigger_pair_weights()`).", call. = FALSE)
  }

  n <- length(bgprobs)
  if (nrow(pairs) > 0L &&
      (any(pairs$i < 1L | pairs$i > n) || any(pairs$j < 1L | pairs$j > n))) {
    stop("`pairs$i` and `pairs$j` must be valid event indices in 1..length(bgprobs).",
         call. = FALSE)
  }
  mark <- rep(0L, n)

  pairs_by_i <- split(seq_len(nrow(pairs)), factor(pairs$i, levels = seq_len(n)))
  for (i in seq_len(n)) {
    u <- stats::runif(1)

    if (u < bgprobs[i]) {
      mark[i] <- 0L
      next
    }

    rows <- pairs_by_i[[i]]

    if (length(rows) > 0) {
      parent_events <- pairs$j[rows]
      rho_values <- pair_weights[rows]

      cumsum_prob <- bgprobs[i] + cumsum(rho_values)
      parent_pos <- which(u <= cumsum_prob)[1]

      mark[i] <- parent_events[if (is.na(parent_pos)) length(parent_events) else parent_pos]
    }
  }

  out <- list(mark = mark)
  class(out) <- "semihawkes_branching_sample"
  out
}

#' @export
print.semihawkes_branching_sample <- function(x, ...) {
  n <- length(x$mark)
  n_bg <- sum(x$mark == 0)
  cat("<semihawkes_branching_sample>\n")
  cat(sprintf("  %d events, %d background (%.1f%%), %d triggered\n",
              n, n_bg, 100 * n_bg / n, n - n_bg))
  invisible(x)
}

#' Draw multiple independent realizations of the branching structure
#'
#' Repeatedly calls `sample_branching()` to build a Monte Carlo ensemble of
#' branching-structure realizations, for quantifying uncertainty in any
#' summary that depends on the background/triggered classification of
#' events (e.g. the number of background events per month, or their
#' spatial concentration -- see the package vignette on branching-structure
#' UQ).
#'
#' @param bgprobs,pairs,pair_weights Passed to `sample_branching()` on
#'   every draw (fixed across draws -- only the internal `u ~ Uniform(0,1)`
#'   draws differ between samples).
#' @param n_samples Number of independent realizations to draw; default
#'   `100`.
#' @param seed Optional seed passed to `set.seed()` before drawing, for
#'   reproducibility. Default `NULL` (no seed set -- uses the current RNG
#'   state).
#' @param restore_rng Only used when `seed` is supplied. If `TRUE` (the
#'   default), the global random-number state is put back as it was before
#'   the call, so the seed affects only this call. If `FALSE`, the call leaves
#'   the generator in the state `set.seed(seed)` and its own draws produce,
#'   as earlier versions of the package did -- useful for reproducing results
#'   from scripts that relied on that.
#'
#' @return A list of length `n_samples`, each element a
#'   `"semihawkes_branching_sample"` object from `sample_branching()`.
#'
#' @seealso [sample_branching()] for a single draw,
#'   [summarize_branching_samples()] to extract per-sample summaries (e.g.
#'   number of background events) from the ensemble.
#'
#' @examples
#' bgprobs <- c(1, 0.4, 0.5, 0.2)
#' pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
#'                     j = as.integer(c(1, 1, 2, 2, 3)))
#' pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)
#'
#' samples <- sample_branching_multi(bgprobs, pairs, pair_weights,
#'                                   n_samples = 200, seed = 1)
#' length(samples)
#' summarize_branching_samples(samples)
#' @export
sample_branching_multi <- function(bgprobs, pairs, pair_weights,
                                   n_samples = 100, seed = NULL,
                                   restore_rng = TRUE) {
  .check_flag(restore_rng, "restore_rng")
  if (!is.null(seed)) {
    if (restore_rng) {
      old_rng <- .get_rng_state()
      on.exit(.restore_rng_state(old_rng), add = TRUE)
    }
    set.seed(seed)
  }

  samples <- vector("list", n_samples)
  for (s in seq_len(n_samples)) {
    samples[[s]] <- sample_branching(bgprobs, pairs, pair_weights)
  }
  samples
}

#' Check consistency of background probabilities and excitation weights
#'
#' For each event `i`, verifies that `bgprobs[i] + sum(pair_weights[pairs$i
#' == i])` equals `1` (up to a numerical tolerance) -- the condition
#' `sample_branching()` relies on for `u ~ Uniform(0, 1)` to be partitioned
#' correctly between the background outcome and every candidate-parent
#' outcome. Call this once before drawing samples with
#' `sample_branching()`/`sample_branching_multi()`.
#'
#' @param bgprobs,pairs,pair_weights Same as in `sample_branching()`.
#' @param tol Numerical tolerance for the sum-to-1 check; default `1e-6`.
#'
#' @return A list of class `"semihawkes_branching_check"` with elements:
#'   \describe{
#'     \item{totals}{Numeric vector, one entry per event: `bgprobs[i] +
#'       sum(pair_weights[pairs$i == i])`.}
#'     \item{max_abs_error}{`max(abs(totals - 1))`.}
#'     \item{n_no_parents}{Number of events with `bgprobs[i] < 1` but no
#'       row in `pairs` (candidate parents) -- these silently fall back to
#'       background in `sample_branching()`, which is only correct if
#'       `bgprobs[i]` is actually `1` for those events.}
#'     \item{ok}{`TRUE` if `max_abs_error <= tol` and `n_no_parents == 0`.}
#'   }
#'
#' @examples
#' bgprobs <- c(1, 0.4, 0.5, 0.2)
#' pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
#'                     j = as.integer(c(1, 1, 2, 2, 3)))
#' pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)
#'
#' check <- check_branching_weights(bgprobs, pairs, pair_weights)
#' check
#' stopifnot(check$ok)
#' @export
check_branching_weights <- function(bgprobs, pairs, pair_weights, tol = 1e-6) {
  n <- length(bgprobs)
  pairs_by_i <- split(seq_len(nrow(pairs)), factor(pairs$i, levels = seq_len(n)))

  totals <- vapply(seq_len(n), function(i) {
    rows <- pairs_by_i[[i]]
    bgprobs[i] + if (length(rows) > 0) sum(pair_weights[rows]) else 0
  }, numeric(1))

  n_no_parents <- sum(vapply(seq_len(n), function(i) {
    length(pairs_by_i[[i]]) == 0 && bgprobs[i] < 1 - tol
  }, logical(1)))

  max_abs_error <- max(abs(totals - 1))

  out <- list(
    totals = totals,
    max_abs_error = max_abs_error,
    n_no_parents = n_no_parents,
    ok = (max_abs_error <= tol) && (n_no_parents == 0)
  )
  class(out) <- "semihawkes_branching_check"
  out
}

#' @export
print.semihawkes_branching_check <- function(x, ...) {
  cat("<semihawkes_branching_check>\n")
  cat(sprintf("  max |total - 1| = %.3g\n", x$max_abs_error))
  cat(sprintf("  events with no candidate parents (bgprobs < 1): %d\n",
              x$n_no_parents))
  cat(sprintf("  ok: %s\n", x$ok))
  invisible(x)
}

#' Summarize an ensemble of branching-structure samples
#'
#' Extracts the number of background events from each realization in a
#' `sample_branching_multi()` ensemble, for reporting the central tendency
#' of the background-event count under sampling uncertainty. Reports the
#' median only (not the mean or any fixed-threshold count such as
#' `sum(bgprobs >= 0.5)`), since the empirical distribution of background
#' counts is often asymmetric and the median is more robust to that
#' asymmetry.
#'
#' @param samples A list of `"semihawkes_branching_sample"` objects, from
#'   `sample_branching_multi()`.
#'
#' @return A list of class `"semihawkes_branching_summary"` with elements:
#'   \describe{
#'     \item{num_bg}{Integer vector, length `length(samples)`: the number
#'       of background events in each sample.}
#'     \item{median_bg}{`median(num_bg)`.}
#'     \item{iqr_bg}{Named numeric vector, `quantile(num_bg, c(0.25,
#'       0.75))`.}
#'   }
#'
#' @examples
#' bgprobs <- c(1, 0.4, 0.5, 0.2)
#' pairs <- data.frame(i = as.integer(c(2, 3, 3, 4, 4)),
#'                     j = as.integer(c(1, 1, 2, 2, 3)))
#' pair_weights <- c(0.6, 0.3, 0.2, 0.5, 0.3)
#'
#' samples <- sample_branching_multi(bgprobs, pairs, pair_weights,
#'                                   n_samples = 200, seed = 1)
#' summarize_branching_samples(samples)
#' @export
summarize_branching_samples <- function(samples) {
  num_bg <- vapply(samples, function(s) sum(s$mark == 0L), integer(1))

  out <- list(
    num_bg = num_bg,
    median_bg = stats::median(num_bg),
    iqr_bg = stats::quantile(num_bg, c(0.25, 0.75))
  )
  class(out) <- "semihawkes_branching_summary"
  out
}

#' @export
print.semihawkes_branching_summary <- function(x, ...) {
  cat("<semihawkes_branching_summary>\n")
  cat(sprintf("  %d samples\n", length(x$num_bg)))
  cat(sprintf("  median background events: %s\n", x$median_bg))
  cat(sprintf("  IQR: [%s, %s]\n", x$iqr_bg[1], x$iqr_bg[2]))
  invisible(x)
}

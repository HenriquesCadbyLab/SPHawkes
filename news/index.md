# Changelog

## spHawkes 0.1.1

First release.

spHawkes fits semi-parametric spatio-temporal Hawkes (self-exciting)
point process models with the stochastic-declustering EM algorithm of
Zhuang and Mateu (2019), on any projected `sf` point layer and
study-region polygon. This release provides:

- data preparation and checks for event points and study regions, with
  adaptive kernel bandwidths and edge-correction integrals
  ([`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md),
  [`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md),
  [`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md)),
  and exploratory plots before fitting
  ([`plot_event_summary()`](https://henriquescadbylab.github.io/spHawkes/reference/plot_event_summary.md),
  [`inspect_spatial_grid()`](https://henriquescadbylab.github.io/spHawkes/reference/inspect_spatial_grid.md));
- the model fit,
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md),
  with a non-parametric background (long-term trend, weekly, monthly and
  spatial components) and non-parametric triggering kernels in time and
  space. Each step of the EM algorithm is also exported on its own, and
  the slowest steps can run on parallel workers (`nprocs`);
- goodness-of-fit diagnostics:
  [`gof_time_transform()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_time_transform.md),
  [`gof_spatial_residuals()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_spatial_residuals.md)
  and
  [`gof_period_residuals()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_period_residuals.md);
- gridded intensity maps over time:
  [`background_trigger_grid()`](https://henriquescadbylab.github.io/spHawkes/reference/background_trigger_grid.md),
  [`plot_intensity_snapshot()`](https://henriquescadbylab.github.io/spHawkes/reference/plot_intensity_snapshot.md),
  [`animate_intensity_snapshots()`](https://henriquescadbylab.github.io/spHawkes/reference/animate_intensity_snapshots.md)
  and
  [`integrate_intensity_period()`](https://henriquescadbylab.github.io/spHawkes/reference/integrate_intensity_period.md);
- Monte Carlo sampling of the branching structure, to quantify how
  events divide between background and triggered:
  [`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md)
  and
  [`summarize_branching_samples()`](https://henriquescadbylab.github.io/spHawkes/reference/summarize_branching_samples.md);
- three example datasets, each with its study region: `sim_events`
  (simulated), `crimes` (simulated from a model fitted to open homicide
  data from Guayaquil, Ecuador) and `milwaukee` (drug-related deaths in
  Milwaukee County, Wisconsin), used in four vignettes.

spHawkes began as a port of the R code that accompanies D’Angelo et
al. (2024). It has since been substantially rewritten: several errors in
that code have been fixed, and its slowest computations have been
rewritten to run much faster. The model and the estimation method are
the same.

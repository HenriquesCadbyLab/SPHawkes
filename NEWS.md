# spHawkes 0.1.1

First release.

spHawkes fits semi-parametric spatio-temporal Hawkes (self-exciting) point process models with the stochastic-declustering EM algorithm of Zhuang and Mateu (2019), on any projected `sf` point layer and study-region polygon. This release provides:

-   data preparation and checks for event points and study regions, with adaptive kernel bandwidths and edge-correction integrals (`prepare_events()`, `add_adaptive_bandwidth()`, `add_bg_integral()`), and exploratory plots before fitting (`plot_event_summary()`, `inspect_spatial_grid()`);
-   the model fit, `fit_semihawkes()`, with a non-parametric background (long-term trend, weekly, monthly and spatial components) and non-parametric triggering kernels in time and space. Each step of the EM algorithm is also exported on its own, and the slowest steps can run on parallel workers (`nprocs`);
-   goodness-of-fit diagnostics: `gof_time_transform()`, `gof_spatial_residuals()` and `gof_period_residuals()`;
-   gridded intensity maps over time: `background_trigger_grid()`, `plot_intensity_snapshot()`, `animate_intensity_snapshots()` and `integrate_intensity_period()`;
-   Monte Carlo sampling of the branching structure, to quantify how events divide between background and triggered: `sample_branching_multi()` and `summarize_branching_samples()`;
-   three example datasets, each with its study region: `sim_events` (simulated), `crimes` (simulated from a model fitted to open homicide data from Guayaquil, Ecuador) and `milwaukee` (drug-related deaths in Milwaukee County, Wisconsin), used in four vignettes.

spHawkes began as a port of the R code that accompanies D'Angelo et al. (2024). It has since been substantially rewritten: several errors in that code have been fixed, and its slowest computations have been rewritten to run much faster. The model and the estimation method are the same.

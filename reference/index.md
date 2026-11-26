# Package index

## Fitting a model

The usual entry point.
[`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
runs the whole pipeline from raw `sf` points to a converged fit; the
stages it composes are exported individually below, so a single step can
be inspected without running the full EM loop.

- [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  : Fit a semi-parametric self-exciting (Hawkes) point process by EM

- [`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md)
  : Single EM step: E-step responsibilities + M-step update of A, mu,
  and all non-parametric background/triggering components

- [`summary(`*`<semihawkes_fit>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/summary.semihawkes_fit.md)
  [`print(`*`<summary.semihawkes_fit>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/summary.semihawkes_fit.md)
  : Summarize a fitted semi-parametric Hawkes model

- [`plot(`*`<semihawkes_fit>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/plot.semihawkes_fit.md)
  :

  Plot method for a fitted `semihawkes_fit` object

## Data preparation

Validating and reshaping event data, and the two per-event quantities
the background estimator needs: an adaptive kernel bandwidth and a Monte
Carlo edge-correction integral.

- [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)
  : Prepare a spatio-temporal event dataset for semi-parametric Hawkes
  fitting

- [`adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/adaptive_bandwidth.md)
  : Adaptive (nearest-neighbour) spatial bandwidth for each event

- [`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md)
  :

  Add an adaptive bandwidth column to a `semihawkes_events` object

- [`bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/bg_integral.md)
  : Monte Carlo integral of each event's isotropic Gaussian kernel over
  the study region

- [`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md)
  :

  Add the background-kernel edge-correction integral to a
  `semihawkes_events` object

## Background intensity

The non-parametric components of the background term — a long-term
trend, calendar-month and day-of-week periodicities, and a spatial
surface.

- [`temporal_base_trend()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_trend.md)
  : Long-term trend temporal base value

- [`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md)
  : Weekly non-parametric temporal base value

- [`temporal_base_monthly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_monthly.md)
  : Monthly non-parametric temporal base value

- [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md)
  : Non-parametric spatial background surface

- [`background_intensity()`](https://henriquescadbylab.github.io/spHawkes/reference/background_intensity.md)
  :

  Combine temporal and spatial background components into the background
  intensity shape (without the scalar `mu` parameter)

## Triggering (self-exciting) component

Initial parametric guesses for the temporal and spatial triggering
kernels, their evaluation over events and the study domain, and their
non-parametric re-estimation from EM responsibility weights.

- [`trigger_temporal_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_temporal_init.md)
  : Initial (parametric) temporal triggering-function guess
- [`trigger_spatial_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_init.md)
  : Initial (parametric) spatial triggering-function guess
- [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md)
  : Evaluate self-exciting (triggering) contributions at events and over
  the study domain
- [`trigger_pairs()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pairs.md)
  : Build the (i, j) triggering-pair index for re-estimation
- [`trigger_pair_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_pair_weights.md)
  : EM responsibility weights for triggering pairs
- [`trigger_spatial_marks()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_marks.md)
  : Cache per-event spatial edge-correction marks
- [`trigger_edge_correction()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_edge_correction.md)
  : Per-event temporal/spatial edge-correction factors for trigger
  re-estimation
- [`trigger_reestimate_temporal()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_temporal.md)
  : Re-estimate the temporal triggering kernel from EM responsibility
  weights
- [`trigger_reestimate_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_spatial.md)
  : Re-estimate the spatial triggering kernel from EM responsibility
  weights

## Branching structure and uncertainty quantification

Stochastic declustering: Monte Carlo sampling of the branching structure
(which events are background vs. triggered, and by whom) from a fitted
model, for quantifying uncertainty in any downstream summary.

- [`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md)
  : Sample one realization of the branching structure via UQ
- [`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md)
  : Draw multiple independent realizations of the branching structure
- [`check_branching_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/check_branching_weights.md)
  : Check consistency of background probabilities and excitation weights
- [`summarize_branching_samples()`](https://henriquescadbylab.github.io/spHawkes/reference/summarize_branching_samples.md)
  : Summarize an ensemble of branching-structure samples

## Diagnostics and visualization

Grid-based intensity snapshots, spatial maps, and animations built from
a fitted model.

- [`plot_event_summary()`](https://henriquescadbylab.github.io/spHawkes/reference/plot_event_summary.md)
  : Plot a summary of a spatio-temporal event data set
- [`background_trigger_grid()`](https://henriquescadbylab.github.io/spHawkes/reference/background_trigger_grid.md)
  : Evaluate background and triggering intensity on a spatial grid over
  time
- [`integrate_intensity_period()`](https://henriquescadbylab.github.io/spHawkes/reference/integrate_intensity_period.md)
  : Integrate gridded intensity over one time period
- [`background_to_grid_sf()`](https://henriquescadbylab.github.io/spHawkes/reference/background_to_grid_sf.md)
  : Convert a semihawkes background grid to an sf polygon grid
- [`crop_background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/crop_background_spatial.md)
  : Crop a semihawkes spatial background to its study-area sf polygon
- [`plot_intensity_snapshot()`](https://henriquescadbylab.github.io/spHawkes/reference/plot_intensity_snapshot.md)
  : Plot one gridded intensity snapshot
- [`animate_intensity_snapshots()`](https://henriquescadbylab.github.io/spHawkes/reference/animate_intensity_snapshots.md)
  : Animate previously evaluated intensity snapshots
- [`inspect_spatial_grid()`](https://henriquescadbylab.github.io/spHawkes/reference/inspect_spatial_grid.md)
  : Inspect a candidate spatial grid

## Simulated example data

A small simulated dataset with known generating values, so that fitted
components have real structure to recover and something to be checked
against.

- [`sim_events`](https://henriquescadbylab.github.io/spHawkes/reference/sim_events.md)
  : Simulated crime events for an example city
- [`sim_region`](https://henriquescadbylab.github.io/spHawkes/reference/sim_region.md)
  : Study-area boundary for the example city

## Milwaukee example data

Geocoded drug-related deaths in Milwaukee County from 2020 onward,
supplied with a polygonal observation window for illustrating
semi-parametric spatio-temporal Hawkes modelling.

- [`milwaukee`](https://henriquescadbylab.github.io/spHawkes/reference/milwaukee.md)
  : Drug-related deaths in Milwaukee County
- [`milwaukee_county`](https://henriquescadbylab.github.io/spHawkes/reference/milwaukee_county.md)
  : Milwaukee County observation window

## Simulated crime example data

Events simulated from a semi-parametric Hawkes model fitted to
transformed intentional homicide records from Guayaquil, Ecuador, with a
geomasked observation window, for illustrating semi-parametric
spatio-temporal Hawkes modelling.

- [`crimes`](https://henriquescadbylab.github.io/spHawkes/reference/crimes.md)
  : Synthetic crime-like events from a fitted Hawkes model
- [`crimes_region`](https://henriquescadbylab.github.io/spHawkes/reference/crimes_region.md)
  : Geomasked spatial window for the crime example

## Numerical primitives

Simpson’s-rule integration, weighted histograms and Gaussian kernel
smoothing, ported from Zhuang & Mateu (2019). Note that the smoothers
differ in their boundary handling and are not interchangeable — see
[`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md).

- [`simpson()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.md)
  : Simpson's rule integral of an equally-spaced function
- [`simpson.2D()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.2D.md)
  : Simpson's rule double integral of an equally-spaced function
- [`hist(`*`<weighted>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md)
  [`hist(`*`<weighted.1D>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.md)
  : Weighted 1D histogram
- [`hist(`*`<weighted.2D>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/hist.weighted.2D.md)
  : Weighted 2D histogram
- [`ker.smooth.conv()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.conv.md)
  : 1D Gaussian kernel smoothing by direct convolution, with reflected
  boundaries
- [`ker.smooth.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md)
  [`ker.smooth.1D.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.fft.md)
  : 1D Gaussian kernel smoothing by FFT, with periodic boundaries
- [`ker.smooth.2D.fft()`](https://henriquescadbylab.github.io/spHawkes/reference/ker.smooth.2D.fft.md)
  : 2D Gaussian kernel smoothing by FFT, with reflected boundaries

## Goodness-of-fit diagnostics

Tools for assessing fitted semi-parametric Hawkes models using temporal
and spatial residual diagnostics.

- [`gof_time_transform()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_time_transform.md)
  : Time-transformed goodness-of-fit diagnostic
- [`gof_spatial_residuals()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_spatial_residuals.md)
  : Spatial raw residuals for a fitted model
- [`gof_period_residuals()`](https://henriquescadbylab.github.io/spHawkes/reference/gof_period_residuals.md)
  : Spatial residuals by time period
- [`intensity_at_times()`](https://henriquescadbylab.github.io/spHawkes/reference/intensity_at_times.md)
  : Total model intensity integrated over space
- [`plot(`*`<semihawkes_gof_time>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/plot.semihawkes_gof_time.md)
  : Plot a time-transformed goodness-of-fit diagnostic
- [`plot(`*`<semihawkes_gof_spatial>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/plot.semihawkes_gof_spatial.md)
  : Plot spatial residuals
- [`plot(`*`<semihawkes_gof_period>`*`)`](https://henriquescadbylab.github.io/spHawkes/reference/plot.semihawkes_gof_period.md)
  : Plot period-wise spatial residuals

## Package overview

- [`spHawkes`](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md)
  [`spHawkes-package`](https://henriquescadbylab.github.io/spHawkes/reference/spHawkes-package.md)
  : spHawkes: Semi-Parametric Self-Exciting Point Process Models

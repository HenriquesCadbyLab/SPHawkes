# spHawkes: Semi-Parametric Self-Exciting Point Process Models

Fits semi-parametric spatio-temporal Hawkes (self-exciting) point
process models to event data via a stochastic-declustering EM algorithm
(Zhuang & Mateu, 2019), on any `sf` point layer and study-region
polygon.

The conditional intensity is the sum of a **background** term, capturing
baseline place-based risk, and a **triggering** term, capturing the
extra risk each past event induces nearby in space and soon after in
time: \$\$\lambda(x, y, t \mid H_t) = \mu_0 \\ \mu_T(t) \\ \mu_m(t) \\
\mu_w(t) \\ \mu_b(x, y) + A \int\_{-\infty}^{t} \int\_{W} g(t - s) \\
h(x - u, y - v) \\ N(du \\ dv \\ ds),\$\$ where \\W\\ is the study
region and \\N\\ is the counting measure of the events, so the
triggering term adds up \\A \\ g(t - t_i) \\ h(x - x_i, y - y_i)\\ over
the events \\i\\ before \\t\\. The background factors (long-term trend,
calendar month, day of week, and a spatial surface) and both triggering
kernels are estimated non-parametrically; \\\mu_0\\ and the branching
ratio \\A\\ are updated in closed form.

## Getting started

[`vignette("getting-started", package = "spHawkes")`](https://henriquescadbylab.github.io/spHawkes/articles/getting-started.md)
walks through the whole workflow on a small example dataset.

The usual entry point is
[`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md),
which runs the full pipeline from raw `sf` points to a converged fit.
The stages it composes are also exported individually, so a single step
can be inspected without running the whole EM loop:

- Data preparation:

  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md),
  [`add_adaptive_bandwidth()`](https://henriquescadbylab.github.io/spHawkes/reference/add_adaptive_bandwidth.md),
  [`add_bg_integral()`](https://henriquescadbylab.github.io/spHawkes/reference/add_bg_integral.md)

- Background components:

  [`temporal_base_trend()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_trend.md),
  [`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md),
  [`temporal_base_monthly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_monthly.md),
  [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md),
  [`background_intensity()`](https://henriquescadbylab.github.io/spHawkes/reference/background_intensity.md)

- Triggering components:

  [`trigger_temporal_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_temporal_init.md),
  [`trigger_spatial_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_init.md),
  [`trigger_evaluate()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_evaluate.md),
  [`trigger_reestimate_temporal()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_temporal.md),
  [`trigger_reestimate_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_reestimate_spatial.md)

- EM:

  [`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md)
  for one iteration,
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  to convergence

- Branching structure & uncertainty:

  [`sample_branching()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching.md),
  [`sample_branching_multi()`](https://henriquescadbylab.github.io/spHawkes/reference/sample_branching_multi.md),
  [`check_branching_weights()`](https://henriquescadbylab.github.io/spHawkes/reference/check_branching_weights.md),
  [`summarize_branching_samples()`](https://henriquescadbylab.github.io/spHawkes/reference/summarize_branching_samples.md)
  – stochastic declustering / Monte Carlo uncertainty quantification for
  which events are background vs. triggered; see
  [`vignette("branching-structure-uq", package = "spHawkes")`](https://henriquescadbylab.github.io/spHawkes/articles/branching-structure-uq.md)

- Diagnostics:

  [`background_trigger_grid()`](https://henriquescadbylab.github.io/spHawkes/reference/background_trigger_grid.md)
  for gridded background/triggering intensity snapshots over time, the
  building block for exceedance-set and persistence analyses

[sim_events](https://henriquescadbylab.github.io/spHawkes/reference/sim_events.md)
and
[sim_region](https://henriquescadbylab.github.io/spHawkes/reference/sim_region.md)
are a small simulated dataset for trying the above out.

## In the code

spHawkes began as a port of the R code that accompanies D'Angelo et al.
(2024), which implements the method of Zhuang and Mateu (2019). It has
since been substantially rewritten: several errors in that code have
been fixed, and its slowest computations have been rewritten to run much
faster. The model and the estimation method are the same.

## Units and scale

Every distance and bandwidth in the package is in the linear units of
the working CRS, and the defaults were tuned for a metre-based
projection of a city-scale dataset.
[`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)
rejects geographic (longitude and latitude) input for this reason, and
warns when the CRS units are not metres. Applying the package to a
different study area means re-choosing `cell_size`, the
triggering-kernel grid extents, and the smoothing bandwidths to suit it;
see the vignette's final section.

## References

Zhuang, J. and Mateu, J. (2019). A semiparametric spatiotemporal
Hawkes-type point process model with periodic background for crime data.
*Journal of the Royal Statistical Society Series A*, 182(3), 919-942.
[doi:10.1111/rssa.12429](https://doi.org/10.1111/rssa.12429)

Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002). Stochastic declustering
of space-time earthquake occurrences. *Journal of the American
Statistical Association*, 97(458), 369-380.
[doi:10.1198/016214502760046925](https://doi.org/10.1198/016214502760046925)

D'Angelo, N., Payares, D., Adelfio, G., and Mateu, J. (2024).
Self-exciting point process modelling of crimes on linear networks.
Statistical Modelling, 24(2), 139–168.
[doi:10.1177/1471082X221094146](https://doi.org/10.1177/1471082X221094146)
.

## See also

Useful links:

- <https://github.com/HenriquesCadbyLab/spHawkes>

- <https://henriquescadbylab.github.io/spHawkes/>

- Report bugs at <https://github.com/HenriquesCadbyLab/spHawkes/issues>

## Author

**Maintainer**: Luis Sotomayor <sotomayorluis441@gmail.com>
([ORCID](https://orcid.org/0009-0003-3903-7388)) \[copyright holder\]

Authors:

- Olatunji Johnson ([ORCID](https://orcid.org/0000-0002-4080-0999))
  \[copyright holder\]

- Ines Henriques-Cadby ([ORCID](https://orcid.org/0000-0003-3916-7556))
  \[copyright holder\]

Other contributors:

- David Buil-Gil \[contributor\]

- Jorge Mateu \[contributor\]

#' spHawkes: Semi-Parametric Self-Exciting Point Process Models
#'
#' @description
#' Fits semi-parametric spatio-temporal Hawkes (self-exciting) point process
#' models to event data via a stochastic-declustering EM algorithm
#' (Zhuang & Mateu, 2019), on any `sf` point layer and study-region polygon.
#'
#' The conditional intensity is the sum of a **background** term, capturing
#' baseline place-based risk, and a **triggering** term, capturing the extra
#' risk each past event induces nearby in space and soon after in time:
#' \deqn{\lambda(x, y, t \mid H_t) =
#'   \mu_0 \, \mu_T(t) \, \mu_m(t) \, \mu_w(t) \, \mu_b(x, y)
#'   + A \int_{-\infty}^{t} \int_{W} g(t - s) \, h(x - u, y - v) \,
#'   N(du \, dv \, ds),}
#' where \eqn{W} is the study region and \eqn{N} is the counting measure of
#' the events, so the triggering term adds up
#' \eqn{A \, g(t - t_i) \, h(x - x_i, y - y_i)} over the events \eqn{i} before
#' \eqn{t}.
#' The background factors (long-term trend, calendar month, day of week, and a
#' spatial surface) and both triggering kernels are estimated
#' non-parametrically; \eqn{\mu_0} and the branching ratio \eqn{A} are updated
#' in closed form.
#'
#' @section Getting started:
#' `vignette("getting-started", package = "spHawkes")` walks through the
#' whole workflow on a small example dataset.
#'
#' The usual entry point is [fit_semihawkes()], which runs the full pipeline
#' from raw `sf` points to a converged fit. The stages it composes are also
#' exported individually, so a single step can be inspected without running
#' the whole EM loop:
#'
#' \describe{
#'   \item{Data preparation}{[prepare_events()],
#'     [add_adaptive_bandwidth()], [add_bg_integral()]}
#'   \item{Background components}{[temporal_base_trend()],
#'     [temporal_base_weekly()], [temporal_base_monthly()],
#'     [background_spatial()], [background_intensity()]}
#'   \item{Triggering components}{[trigger_temporal_init()],
#'     [trigger_spatial_init()], [trigger_evaluate()],
#'     [trigger_reestimate_temporal()], [trigger_reestimate_spatial()]}
#'   \item{EM}{[em_step()] for one iteration, [fit_semihawkes()] to
#'     convergence}
#'   \item{Branching structure & uncertainty}{[sample_branching()],
#'     [sample_branching_multi()], [check_branching_weights()],
#'     [summarize_branching_samples()] -- stochastic declustering /
#'     Monte Carlo uncertainty quantification for which events are
#'     background vs. triggered; see
#'     `vignette("branching-structure-uq", package = "spHawkes")`}
#'   \item{Diagnostics}{[background_trigger_grid()] for gridded
#'     background/triggering intensity snapshots over time, the building
#'     block for exceedance-set and persistence analyses}
#' }
#'
#' [sim_events] and [sim_region] are a small simulated dataset for trying
#' the above out.
#'
#' @section In the code:
#' spHawkes began as a port of the R code that accompanies D'Angelo et al.
#' (2024), which implements the method of Zhuang and Mateu (2019). It has
#' since been substantially rewritten: several errors in that code have been
#' fixed, and its slowest computations have been rewritten to run much
#' faster. The model and the estimation method are the same.
#'
#' @section Units and scale:
#' Every distance and bandwidth in the package is in the linear units of the
#' working CRS, and the defaults were tuned for a metre-based projection of a
#' city-scale dataset. [prepare_events()] rejects geographic (longitude and
#' latitude) input for this reason, and warns when the CRS units are not
#' metres. Applying the package to a different study area means re-choosing
#' `cell_size`, the triggering-kernel grid extents, and the smoothing
#' bandwidths to suit it; see the vignette's final section.
#'
#' @references
#' Zhuang, J. and Mateu, J. (2019). A semiparametric spatiotemporal
#' Hawkes-type point process model with periodic background for crime data.
#' *Journal of the Royal Statistical Society Series A*, 182(3), 919-942.
#' \doi{10.1111/rssa.12429}
#'
#' Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002). Stochastic declustering of
#' space-time earthquake occurrences. *Journal of the American Statistical
#' Association*, 97(458), 369-380. \doi{10.1198/016214502760046925}
#'
#' D'Angelo, N., Payares, D., Adelfio, G., and Mateu, J. (2024).
#' Self-exciting point process modelling of crimes on linear networks.
#' Statistical Modelling, 24(2), 139--168.
#' \doi{10.1177/1471082X221094146}.
#'
#' @keywords internal
"_PACKAGE"

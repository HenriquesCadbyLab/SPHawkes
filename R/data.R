#' Simulated crime events for an example city
#'
#' A small simulated spatio-temporal point pattern, provided so that the
#' package's examples and vignette can be run without a real crime dataset.
#'
#' **These are not real crime records.** The pattern is drawn from a branching
#' process shaped like the model [fit_semihawkes()] estimates, so each fitted
#' component has genuine structure to recover rather than noise:
#'
#' \describe{
#'   \item{spatial background}{three Gaussian hotspots over a uniform floor}
#'   \item{long-term trend}{immigration rate rising across the three years}
#'   \item{day of week}{weekends elevated relative to weekdays}
#'   \item{triggering}{each event spawns `Poisson(0.35)` offspring, roughly
#'     300 m and 5 days away}
#' }
#'
#' The generating values are recorded here because they are useful as targets
#' when checking that the estimator recovers what it should: branching ratio
#' `A = 0.35`, triggering length scale 300 m, triggering time scale 5 days.
#' Of the 375 simulated events, 260 are background and 115 (30.7%) triggered.
#'
#' @format An [sf::sf] object with 375 rows and 3 columns:
#' \describe{
#'   \item{id}{Integer event identifier.}
#'   \item{date}{Event date, a `Date`, between 2021-01-01 and 2023-12-31.}
#'   \item{geometry}{`POINT` geometry in EPSG:32717 (UTM zone 17S, metres).}
#' }
#'
#' @source Simulated by `data-raw/make-example-data.R` in the package's GitHub
#'   repository, <https://github.com/HenriquesCadbyLab/spHawkes>.
#' @seealso [sim_region] for the matching study-area boundary.
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' print(ev)
"sim_events"

#' Study-area boundary for the example city
#'
#' The square study region that [sim_events] was simulated within: 8 km by
#' 8 km, in EPSG:32717 (UTM zone 17S), so all coordinates and distances are in
#' metres. This is the projected, metre-based setting the package's default
#' bandwidths and grid extents assume.
#'
#' @format An [sf::sf] object with 1 row and 2 columns:
#' \describe{
#'   \item{name}{Region name, a character string.}
#'   \item{geometry}{`POLYGON` geometry in EPSG:32717.}
#' }
#'
#' @source Simulated by `data-raw/make-example-data.R` in the package's GitHub
#'   repository, <https://github.com/HenriquesCadbyLab/spHawkes>.
#' @seealso [sim_events] for the matching event locations.
"sim_region"

#' Synthetic crime-like events from a fitted Hawkes model
#'
#' Simulated spatial-temporal point events provided to demonstrate
#' \code{\link{fit_semihawkes}}. The events were generated from a
#' semi-parametric Hawkes model fitted to previously transformed crime data
#' from Guayaquil, Ecuador.
#' Each row represents a simulated event, not an individual source record.
#'
#' @details
#' The simulation uses the fitted spatial background, long-term trend,
#' weekly and monthly background components, and spatial and temporal
#' triggering kernels.
#'
#' Background events and their descendants were generated using a branching
#' construction with Poisson thinning. The simulation starts with an empty
#' event history. Descendants falling outside the observation window are
#' discarded; events outside the window do not generate further descendants
#' that can enter it.
#'
#' Coordinates use a local Cartesian reference system, with distances
#' expressed in metres. They do not identify the original geographic
#' locations. Calendar dates belong to the transformed modelling timeline
#' and should not be interpreted as dates of recorded crimes.
#'
#' The dataset is intended for software examples and methodological
#' illustration. It should not be used to draw substantive conclusions
#' about crime patterns in the source population.
#'
#' @format An \code{sf} object with the following variables:
#' \describe{
#'   \item{id}{Integer. Sequential identifier for each simulated event,
#'   ordered by event time.}
#'   \item{date}{Simulated event date-time, stored as \code{POSIXct} in UTC.}
#'   \item{geometry}{Simulated \code{POINT} geometry in local Cartesian
#'   coordinates, expressed in metres.}
#' }
#'
#' @source
#' Generated for \pkg{spHawkes} by simulating from a semi-parametric Hawkes
#' model, as described under Details. The model was fitted to intentional
#' homicide records for Guayaquil, published as open data by Ecuador's
#' Ministry of the Interior, after their coordinates and dates had been
#' transformed. The
#' packaged events are newly simulated rather than copies of those records.
#'
#' Source of the records: Ministerio del Interior del Ecuador,
#' \enc{Subsecretaría}{Subsecretaria} de Estudios y
#' \enc{Estadística}{Estadistica} de la Seguridad. Homicidios Intencionales.
#' Datos Abiertos Ecuador.
#' <https://www.datosabiertos.gob.ec/dataset/homicidios-intencionales?activity_id=0da52f2a-0a1e-4d0d-afa4-0fae1adab14b>.
#' Licensed under Creative Commons Attribution (CC BY).
#'
#' @seealso \code{\link{crimes_region}}, \code{\link{fit_semihawkes}}
#'
#' @examples
#' data("crimes")
#' head(crimes)
#'
#' @name crimes
#' @docType data
#' @keywords datasets
"crimes"


#' Geomasked spatial window for the crime example
#'
#' A modified spatial window accompanying \code{\link{crimes}}, a dataset of
#' events simulated from a model fitted to crime data from Guayaquil,
#' Ecuador.
#'
#' It is provided
#' for methodological demonstration and should not be interpreted
#' as the original administrative boundary of the city. For the source and
#' licence of the records the \code{\link{crimes}} events were simulated
#' from, see \code{\link{crimes}}.
#'
#' @format An \code{sf} object containing a geometry column in a
#'   local Cartesian reference system, with distances expressed
#'   in metres.
#'
#' @seealso \code{\link{crimes}}
#'
#' @examples
#' data("crimes")
#' data("crimes_region")
#'
#' plot(sf::st_geometry(crimes_region), col = "grey95")
#' plot(
#'   sf::st_geometry(crimes),
#'   add = TRUE,
#'   pch = 16,
#'   cex = 0.3
#' )
"crimes_region"


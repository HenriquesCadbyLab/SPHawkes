#' Plot a summary of a spatio-temporal event data set
#'
#' Creates a three-panel diagnostic display containing the event locations,
#' gridded spatial counts, and cumulative event count.
#'
#' @param x A `semihawkes_events` object, as returned by [prepare_events()].
#' @param region An `sf` polygon defining the study region. Defaults to
#'   `x$region`.
#' @param start_date,end_date Study bounds; `end_date` is exclusive. They
#'   default to `x$start_date` and `x$end_date`.
#' @param grid_n Approximate number of square cells across the longest side of
#'   the region's bounding box.
#'
#' @return A `patchwork` object.
#' @importFrom rlang .data
#' @examples
#' ev <- prepare_events(sim_events, sim_region, date_col = "date")
#' plot_event_summary(ev)
#' @export
plot_event_summary <- function(x,
                               region = x$region,
                               start_date = x$start_date,
                               end_date = x$end_date,
                               grid_n = 10L) {
  if (!inherits(x, "semihawkes_events")) {
    stop("`x` must be a `semihawkes_events` object from `prepare_events()`.",
         call. = FALSE)
  }
  if (!inherits(region, c("sf", "sfc"))) {
    stop("`region` must be an sf polygon.", call. = FALSE)
  }
  if (isTRUE(sf::st_is_longlat(region))) {
    stop("A projected CRS is required for square grids.", call. = FALSE)
  }

  dates <- as.Date(x$events$date)
  if (anyNA(dates)) stop("`date` contains missing or invalid values.",
                         call. = FALSE)
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  if (end_date <= start_date || any(dates < start_date | dates >= end_date)) {
    stop("Events must lie in [`start_date`, `end_date`).", call. = FALSE)
  }

  # Region and event locations ---------------------------------------------
  region_sf <- sf::st_sf(
    geometry = sf::st_union(sf::st_make_valid(sf::st_geometry(region)))
  )
  points <- sf::st_as_sf(x$events, coords = c("coorx", "coory"),
                         crs = sf::st_crs(region_sf))
  area_km2 <- as.numeric(units::set_units(sum(sf::st_area(region_sf)), "km^2"))

  # Daily counts and cumulative frequency -----------------------------------
  all_days <- seq(start_date, end_date - 1, by = "day")
  daily <- data.frame(
    date = all_days,
    count = tabulate(match(dates, all_days), nbins = length(all_days))
  )
  daily$cumulative <- cumsum(daily$count)

  # Spatial grid ------------------------------------------------------------
  bbox <- sf::st_bbox(region_sf)
  cellsize <- max(bbox[["xmax"]] - bbox[["xmin"]],
                  bbox[["ymax"]] - bbox[["ymin"]]) / grid_n
  grid <- sf::st_make_grid(region_sf, cellsize = cellsize, square = TRUE)
  grid <- suppressWarnings(sf::st_intersection(sf::st_sf(geometry = grid),
                                               region_sf))
  grid <- grid[sf::st_geometry_type(grid) %in% c("POLYGON", "MULTIPOLYGON"), ]
  grid$count <- lengths(sf::st_intersects(grid, points))

  # Draw in map units. coord_sf() transforms to lon/lat to build the axes,
  # which fails for local (engineering) CRSs, so the CRS is dropped here.
  region_sf <- sf::st_set_crs(region_sf, NA)
  points <- sf::st_set_crs(points, NA)
  grid <- sf::st_set_crs(grid, NA)

  # Panels ------------------------------------------------------------------
  ink <- "#253746"
  outline <- "#626262"

  # Identical spatial limits and padding for both maps
  xlim <- as.numeric(bbox[c("xmin", "xmax")])
  ylim <- as.numeric(bbox[c("ymin", "ymax")])

  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.025
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.025

  # Shared typography and spacing
  theme_common <- ggplot2::theme(
    text = ggplot2::element_text(
      family = "sans",
      size = 9,
      colour = ink
    ),
    plot.title = ggplot2::element_text(
      size = 9,
      face = "bold",
      hjust = 0,
      margin = ggplot2::margin(b = 5)
    ),
    plot.title.position = "plot",
    plot.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

  theme_map <- ggplot2::theme_void(
    base_size = 9,
    base_family = "sans"
  ) +
    theme_common +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(1, "mm")
    )

  p_locations <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = region_sf,
      fill = "#F7F8FA",
      colour = outline,
      linewidth = 0.3
    ) +
    ggplot2::geom_sf(
      data = points,
      colour = ink,
      shape = 16,
      size = 0.4,
      alpha = 0.45
    ) +
    ggplot2::coord_sf(
      datum = NA,
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::labs(title = "(a) Event locations") +
    theme_map

  p_grid <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = grid,
      ggplot2::aes(fill = .data$count),
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::geom_sf(
      data = region_sf,
      fill = NA,
      colour = outline,
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c(
        "#F7FBFF",
        "#DEEBF7",
        "#9ECAE1",
        "#4292C6",
        "#08519C"
      ),
      trans = "sqrt",
      name = "Events",
      breaks = function(z) {
        b <- pretty(z, n = 3)
        b[b >= 0 & b == floor(b)]
      },
      labels = function(z) {
        format(z, big.mark = ",", trim = TRUE, scientific = FALSE)
      },
      guide = ggplot2::guide_colourbar(
        direction = "horizontal",
        theme = ggplot2::theme(
          legend.key.width = grid::unit(45, "mm"),
          legend.key.height = grid::unit(2, "mm"),
          legend.title.position = "top",
          legend.text = ggplot2::element_text(size = 7)
        )
      )
    ) +
    ggplot2::coord_sf(
      datum = NA,
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    ggplot2::labs(title = "(b) Event count") +
    theme_map

  p_cumulative <- ggplot2::ggplot(
    daily,
    ggplot2::aes(x = .data$date, y = .data$cumulative)
  ) +
    ggplot2::geom_step(
      colour = ink,
      linewidth = 0.5
    ) +
    ggplot2::scale_x_date(
      limits = c(start_date, end_date),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      n.breaks = 4,
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      labels = function(z) {
        format(z, big.mark = ",", trim = TRUE, scientific = FALSE)
      }
    ) +
    ggplot2::labs(
      title = "(c) Cumulative events",
      x = NULL,
      y = "Events"
    ) +
    ggplot2::theme_classic(
      base_size = 9,
      base_family = "sans"
    ) +
    theme_common +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        size = 8,
        colour = "#4D4D4D"
      ),
      axis.title.y = ggplot2::element_text(
        size = 8,
        margin = ggplot2::margin(r = 6)
      ),
      axis.line = ggplot2::element_line(
        colour = outline,
        linewidth = 0.25
      ),
      axis.ticks = ggplot2::element_line(
        colour = outline,
        linewidth = 0.25
      ),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.margin = ggplot2::margin(8, 7, 5, 7)
    )

  patchwork::wrap_plots(
    p_locations,
    p_grid,
    p_cumulative,
    design = "AB\nCC",
    heights = c(3, 1),
    guides = "keep"
  ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        plot.margin = ggplot2::margin(3, 3, 3, 3)
      )
    )
}


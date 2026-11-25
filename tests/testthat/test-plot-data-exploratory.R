testthat::test_that("plot_event_summary returns a patchwork without warnings", {
  events <- prepare_events(
    points = sim_events,
    region = sim_region,
    date_col = "date",
    crs = sf::st_crs(sim_region)
  )

  testthat::expect_no_warning({
    p <- plot_event_summary(events)
    patchwork::patchworkGrob(p)
  })

  testthat::expect_s3_class(p, "patchwork")
})

testthat::test_that("plot_event_summary accepts configurable arguments", {
  events <- prepare_events(
    points = sim_events,
    region = sim_region,
    date_col = "date",
    crs = sf::st_crs(sim_region)
  )

  testthat::expect_no_error(
    plot_event_summary(
      events,
      region = sim_region,
      start_date = events$start_date,
      end_date = events$end_date,
      grid_n = 15
    )
  )
})

testthat::test_that("plot_event_summary rejects invalid inputs", {
  events <- prepare_events(
    points = sim_events,
    region = sim_region,
    date_col = "date",
    crs = sf::st_crs(sim_region)
  )

  testthat::expect_error(
    plot_event_summary(sim_events),
    "`x` must be a `semihawkes_events` object",
    fixed = TRUE
  )

  testthat::expect_error(
    plot_event_summary(events, region = data.frame()),
    "`region` must be an sf polygon.",
    fixed = TRUE
  )

  testthat::expect_error(
    plot_event_summary(
      events,
      region = sf::st_transform(sim_region, 4326)
    ),
    "A projected CRS is required",
    fixed = TRUE
  )
})

testthat::test_that("plot_event_summary validates event dates", {
  events <- prepare_events(
    points = sim_events,
    region = sim_region,
    date_col = "date",
    crs = sf::st_crs(sim_region)
  )

  missing_date <- events
  missing_date$events$date[1] <- as.Date(NA)

  testthat::expect_error(
    plot_event_summary(missing_date),
    "`date` contains missing or invalid values.",
    fixed = TRUE
  )

  outside_period <- events
  outside_period$events$date[1] <- outside_period$end_date

  testthat::expect_error(
    plot_event_summary(outside_period),
    "Events must lie in [`start_date`, `end_date`).",
    fixed = TRUE
  )
})

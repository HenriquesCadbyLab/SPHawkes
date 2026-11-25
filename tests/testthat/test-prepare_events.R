test_that("prepare_events validates inputs", {
  d <- make_synthetic_data()

  expect_error(prepare_events(d$points$geometry, d$region, "event_date"),
               "must be an object of class 'sf'")
  expect_error(prepare_events(sf::st_drop_geometry(d$points), d$region, "event_date"),
               "must be an object of class 'sf'")
  expect_error(prepare_events(d$points, d$region, "not_a_column"),
               "not a column")
})

test_that("prepare_events rejects geographic CRS", {
  d <- make_synthetic_data()
  ll_points <- sf::st_transform(d$points, crs = 4326)
  ll_region <- sf::st_transform(d$region, crs = 4326)
  expect_error(prepare_events(ll_points, ll_region, "event_date"),
               "geographic")
})

test_that("prepare_events reprojects points to region CRS", {
  d <- make_synthetic_data()
  ll_points <- sf::st_transform(d$points, crs = 4326) # unprojected points
  out <- prepare_events(ll_points, d$region, "event_date") # region stays projected
  expect_s3_class(out, "semihawkes_events")
  expect_true(all(c("coorx", "coory", "days", "month", "date") %in% names(out$events)))
})

test_that("prepare_events accepts data already in a local engineering CRS", {
  # `crimes` uses a "Local Cartesian coordinates" engineering CRS, for which
  # GDAL/PROJ cannot build a transformation even to itself, so this only
  # works if prepare_events() skips the no-op reprojection.
  out <- prepare_events(crimes, crimes_region, "date")
  expect_s3_class(out, "semihawkes_events")
  expect_equal(out$n_kept, nrow(crimes))
  expect_equal(out$crs, sf::st_crs(crimes_region))

  out2 <- prepare_events(crimes, crimes_region, "date",
                         crs = sf::st_crs(crimes_region))
  expect_equal(out2$n_kept, nrow(crimes))
})

test_that("prepare_events produces events inside the region and date range", {
  d <- make_synthetic_data(n = 300)
  out <- prepare_events(d$points, d$region, "event_date",
                         start_date = "2020-02-01", end_date = "2020-06-01")

  expect_true(all(out$events$days >= 0))
  expect_true(all(out$events$date >= as.Date("2020-02-01")))
  expect_true(all(out$events$date < as.Date("2020-06-01")))
  expect_equal(out$TT, max(out$events$days) + 1)
  expect_true(!is.unsorted(out$events$days))
  expect_equal(out$n_kept, nrow(out$events))
})

test_that("prepare_events converts date-times in their own time zone", {
  # as.Date() on a POSIXct uses UTC, so a Friday-evening event in Chicago
  # (02:00 on Saturday in UTC) used to be recorded on the Saturday.
  synth <- make_synthetic_data(n = 3, seed = 1)
  pts <- synth$points
  pts$when <- as.POSIXct(c("2021-03-05 20:00", "2021-03-06 09:00", "2021-03-08 23:30"),
                         tz = "America/Chicago")
  ev <- prepare_events(pts, synth$region, date_col = "when")
  expect_equal(sort(ev$events$date), as.Date(c("2021-03-05", "2021-03-06", "2021-03-08")))
  expect_equal(ev$start_date, as.Date("2021-03-05"))

  # East of Greenwich the shift went the other way: 01:00 in Tokyo is still
  # the previous day in UTC.
  pts$when <- as.POSIXct(c("2021-03-05 01:00", "2021-03-06 12:00", "2021-03-08 08:30"),
                         tz = "Asia/Tokyo")
  ev <- prepare_events(pts, synth$region, date_col = "when")
  expect_equal(sort(ev$events$date), as.Date(c("2021-03-05", "2021-03-06", "2021-03-08")))

  # Date-time study bounds are converted the same way.
  pts$when <- as.POSIXct(c("2021-03-05 20:00", "2021-03-06 09:00", "2021-03-08 23:30"),
                         tz = "America/Chicago")
  ev <- prepare_events(pts, synth$region, date_col = "when",
                       start_date = as.POSIXct("2021-03-06 20:00", tz = "America/Chicago"))
  expect_equal(sort(ev$events$date), as.Date(c("2021-03-06", "2021-03-08")))

  # UTC date-times, like the bundled `crimes`, give the same dates as before.
  pts$when <- as.POSIXct(c("2021-03-05 20:00", "2021-03-06 09:00", "2021-03-08 23:30"),
                         tz = "UTC")
  ev <- prepare_events(pts, synth$region, date_col = "when")
  expect_equal(sort(ev$events$date), as.Date(sort(pts$when)))
})

test_that("prepare_events errors when nothing survives filtering", {
  d <- make_synthetic_data()
  expect_error(
    prepare_events(d$points, d$region, "event_date",
                    start_date = "2099-01-01", end_date = "2099-02-01"),
    "No events remain"
  )
})

test_that("print/summary methods run without error", {
  d <- make_synthetic_data()
  out <- prepare_events(d$points, d$region, "event_date")
  expect_output(print(out), "semihawkes_events")
  expect_output(summary(out), "semihawkes_events")
})

test_that("the modelled period ends the day after the last event, whatever end_date is", {
  d <- make_synthetic_data(n = 300)
  out3 <- prepare_events(d$points, d$region, "event_date")
  last <- max(out3$events$date)

  # A later `end_date` selects events and labels outputs, but does not move
  # the end of the modelled period; a message says so.
  expect_message(
    out <- prepare_events(d$points, d$region, "event_date",
                          end_date = last + 30),
    "modelled study period ends on"
  )
  expect_equal(out$end_date, last + 30)
  expect_equal(out$start_date + out$TT, last + 1)
  expect_match(paste(capture.output(print(out)), collapse = "\n"),
               paste0("modelled to ", format(last + 1)))

  # No message when the two agree, and the default `end_date` follows the
  # last event kept.
  expect_silent(prepare_events(d$points, d$region, "event_date",
                               end_date = last + 1))
  expect_equal(out3$end_date, out3$start_date + out3$TT)
  expect_no_match(paste(capture.output(print(out3)), collapse = "\n"),
                  "modelled to")
})

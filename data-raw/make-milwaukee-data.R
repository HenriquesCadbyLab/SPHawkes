# data-raw/make-milwaukee-data.R
#
# Construct the Milwaukee County example data distributed with spHawkes.
#
# Source data:
# - Milwaukee County Medical Examiner records and geocoded locations,
#   distributed through the UK-IPOP/open-data-pipeline GitHub release
#   dated 2025-12-01:
#   https://github.com/UK-IPOP/open-data-pipeline/releases/tag/2025-12-01
# - Milwaukee County boundary from U.S. Census Bureau TIGER/Line,
#   obtained with tigris::counties(state = "WI", cb = TRUE, year = 2024).
#
# Construction:
# - Retain Drug Related deaths dated on or after 2020-01-01.
# - Retain records with non-missing geocoded longitude and latitude.
# - Treat source coordinates as WGS84 longitude/latitude (EPSG:4326).
# - Transform events and boundary to NAD83 / UTM zone 16N (EPSG:26916).
# - Retain only events intersecting the Milwaukee County boundary.
#
# The resulting objects are demonstration data for spHawkes.
# They are not an authoritative administrative release.

library(readr)
library(dplyr)
library(sf)
library(tigris)
library(lubridate)

source_url <- "https://github.com/UK-IPOP/open-data-pipeline/releases/download/2025-12-01/Milwaukee_County.zip"

download_dir <- file.path(tempdir(), "spHawkes-milwaukee-data")

dir.create(
  download_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_file <- file.path(
  download_dir,
  basename(source_url)
)

utils::download.file(
  url = source_url,
  destfile = archive_file,
  mode = "wb",
  quiet = FALSE
)

utils::unzip(
  zipfile = archive_file,
  exdir = download_dir
)

records_file <- list.files(
  download_dir,
  pattern = "^records\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

geocoded_file <- list.files(
  download_dir,
  pattern = "^geocoded\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(records_file) != 1L || length(geocoded_file) != 1L) {
  stop(
    "Expected exactly one records.csv and one geocoded.csv after ",
    "extracting the release ZIP.",
    call. = FALSE
  )
}

records <- read_csv(
  records_file,
  show_col_types = FALSE
)

geocoded <- read_csv(
  geocoded_file,
  show_col_types = FALSE
)

records_geocoded <- records %>%
  mutate(
    DeathDate = as.Date(DeathDate)
  ) %>%
  filter(
    !is.na(DeathDate),
    DeathDate >= as.Date("2020-01-01"),
    DeathType == "Drug Related"
  ) %>%
  inner_join(
    geocoded,
    by = "CaseIdentifier"
  ) %>%
  filter(
    !is.na(geocoded_longitude),
    !is.na(geocoded_latitude)
  ) %>%
  select(
    CaseIdentifier,
    DeathDate,
    geocoded_longitude,
    geocoded_latitude
  ) %>%
  st_as_sf(
    coords = c("geocoded_longitude", "geocoded_latitude"),
    crs = 4326,
    remove = TRUE
  ) %>%
  st_transform(crs = 26916)

options(tigris_use_cache = TRUE)

milwaukee_county <- counties(
  state = "WI",
  cb = TRUE,
  year = 2024,
  class = "sf"
) %>%
  filter(NAME == "Milwaukee") %>%
  st_transform(crs = 26916) %>%
  transmute(
    region_name = "Milwaukee County, WI",
    geometry
  )

milwaukee <- records_geocoded[
  lengths(
    st_intersects(records_geocoded, milwaukee_county)
  ) > 0,
] %>%
  mutate(
    event_id = seq_len(n()),
    death_year = year(DeathDate),
    death_month = floor_date(
      DeathDate,
      unit = "month"
    ),
    death_week = floor_date(
      DeathDate,
      unit = "week",
      week_start = 1
    ),
    death_day_of_week = wday(
      DeathDate,
      label = TRUE,
      abbr = TRUE,
      week_start = 1
    ),
    death_weekend = death_day_of_week %in% c("Sat", "Sun")
  ) %>%
  select(
    event_id,
    DeathDate,
    death_year,
    death_month,
    death_week,
    death_day_of_week,
    death_weekend,
    geometry
  )

stopifnot(
  inherits(milwaukee, "sf"),
  inherits(milwaukee_county, "sf"),
  inherits(milwaukee$DeathDate, "Date"),
  all(!is.na(milwaukee$DeathDate)),
  all(milwaukee$DeathDate >= as.Date("2020-01-01")),
  all(st_geometry_type(milwaukee) == "POINT"),
  all(
    st_geometry_type(milwaukee_county) %in% c(
      "POLYGON",
      "MULTIPOLYGON"
    )
  ),
  isTRUE(st_crs(milwaukee)$epsg == 26916),
  isTRUE(st_crs(milwaukee_county)$epsg == 26916),
  all(
    lengths(
      st_intersects(milwaukee, milwaukee_county)
    ) > 0
  )
)

usethis::use_data(
  milwaukee,
  milwaukee_county,
  overwrite = TRUE
)

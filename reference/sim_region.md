# Study-area boundary for the example city

The square study region that
[sim_events](https://henriquescadbylab.github.io/spHawkes/reference/sim_events.md)
was simulated within: 8 km by 8 km, in EPSG:32717 (UTM zone 17S), so all
coordinates and distances are in metres. This is the projected,
metre-based setting the package's default bandwidths and grid extents
assume.

## Usage

``` r
sim_region
```

## Format

An [sf::sf](https://r-spatial.github.io/sf/reference/sf.html) object
with 1 row and 2 columns:

- name:

  Region name, a character string.

- geometry:

  `POLYGON` geometry in EPSG:32717.

## Source

Simulated by `data-raw/make-example-data.R` in the package's GitHub
repository, <https://github.com/HenriquesCadbyLab/spHawkes>.

## See also

[sim_events](https://henriquescadbylab.github.io/spHawkes/reference/sim_events.md)
for the matching event locations.

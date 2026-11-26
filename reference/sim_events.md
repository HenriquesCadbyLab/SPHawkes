# Simulated crime events for an example city

A small simulated spatio-temporal point pattern, provided so that the
package's examples and vignette can be run without a real crime dataset.

## Usage

``` r
sim_events
```

## Format

An [sf::sf](https://r-spatial.github.io/sf/reference/sf.html) object
with 375 rows and 3 columns:

- id:

  Integer event identifier.

- date:

  Event date, a `Date`, between 2021-01-01 and 2023-12-31.

- geometry:

  `POINT` geometry in EPSG:32717 (UTM zone 17S, metres).

## Source

Simulated by `data-raw/make-example-data.R` in the package's GitHub
repository, <https://github.com/HenriquesCadbyLab/spHawkes>.

## Details

**These are not real crime records.** The pattern is drawn from a
branching process shaped like the model
[`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
estimates, so each fitted component has genuine structure to recover
rather than noise:

- spatial background:

  three Gaussian hotspots over a uniform floor

- long-term trend:

  immigration rate rising across the three years

- day of week:

  weekends elevated relative to weekdays

- triggering:

  each event spawns `Poisson(0.35)` offspring, roughly 300 m and 5 days
  away

The generating values are recorded here because they are useful as
targets when checking that the estimator recovers what it should:
branching ratio `A = 0.35`, triggering length scale 300 m, triggering
time scale 5 days. Of the 375 simulated events, 260 are background and
115 (30.7%) triggered.

## See also

[sim_region](https://henriquescadbylab.github.io/spHawkes/reference/sim_region.md)
for the matching study-area boundary.

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date")
print(ev)
#> <semihawkes_events>
#>   Study period : 2021-01-03 to 2024-01-01 (TT = 1093.0 days)
#>   Events kept  : 375 / 375 input rows
#>   Dropped      : 0 outside region, 0 outside date range
#>   CRS          : EPSG:32717
```

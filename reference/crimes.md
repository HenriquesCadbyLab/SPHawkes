# Synthetic crime-like events from a fitted Hawkes model

Simulated spatial-temporal point events provided to demonstrate
[`fit_semihawkes`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).
The events were generated from a semi-parametric Hawkes model fitted to
previously transformed crime data from Guayaquil, Ecuador. Each row
represents a simulated event, not an individual source record.

## Usage

``` r
crimes
```

## Format

An `sf` object with the following variables:

- id:

  Integer. Sequential identifier for each simulated event, ordered by
  event time.

- date:

  Simulated event date-time, stored as `POSIXct` in UTC.

- geometry:

  Simulated `POINT` geometry in local Cartesian coordinates, expressed
  in metres.

## Source

Generated for spHawkes by simulating from a semi-parametric Hawkes
model, as described under Details. The model was fitted to intentional
homicide records for Guayaquil, published as open data by Ecuador's
Ministry of the Interior, after their coordinates and dates had been
transformed. The packaged events are newly simulated rather than copies
of those records.

Source of the records: Ministerio del Interior del Ecuador,
Subsecretaría de Estudios y Estadística de la Seguridad. Homicidios
Intencionales. Datos Abiertos Ecuador.
<https://www.datosabiertos.gob.ec/dataset/homicidios-intencionales?activity_id=0da52f2a-0a1e-4d0d-afa4-0fae1adab14b>.
Licensed under Creative Commons Attribution (CC BY).

## Details

The simulation uses the fitted spatial background, long-term trend,
weekly and monthly background components, and spatial and temporal
triggering kernels.

Background events and their descendants were generated using a branching
construction with Poisson thinning. The simulation starts with an empty
event history. Descendants falling outside the observation window are
discarded; events outside the window do not generate further descendants
that can enter it.

Coordinates use a local Cartesian reference system, with distances
expressed in metres. They do not identify the original geographic
locations. Calendar dates belong to the transformed modelling timeline
and should not be interpreted as dates of recorded crimes.

The dataset is intended for software examples and methodological
illustration. It should not be used to draw substantive conclusions
about crime patterns in the source population.

## See also

[`crimes_region`](https://henriquescadbylab.github.io/spHawkes/reference/crimes_region.md),
[`fit_semihawkes`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)

## Examples

``` r
data("crimes")
head(crimes)
#> Simple feature collection with 6 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 4908.171 ymin: 3239.169 xmax: 9871.017 ymax: 6491.946
#> Projected CRS: Local Cartesian coordinates
#>   id                date                  geometry
#> 1  1 2001-01-07 23:37:17 POINT (6360.513 6169.554)
#> 2  2 2001-01-10 15:49:32 POINT (8133.878 6085.512)
#> 3  3 2001-01-17 07:50:08 POINT (4908.171 3239.169)
#> 4  4 2001-01-19 09:10:29 POINT (6615.126 6491.946)
#> 5  5 2001-01-20 00:58:20 POINT (9433.131 6339.906)
#> 6  6 2001-01-22 18:28:19 POINT (9871.017 5552.274)
```

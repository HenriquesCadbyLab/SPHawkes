# Geomasked spatial window for the crime example

A modified spatial window accompanying
[`crimes`](https://henriquescadbylab.github.io/spHawkes/reference/crimes.md),
a dataset of events simulated from a model fitted to crime data from
Guayaquil, Ecuador.

## Usage

``` r
crimes_region
```

## Format

An `sf` object containing a geometry column in a local Cartesian
reference system, with distances expressed in metres.

## Details

It is provided for methodological demonstration and should not be
interpreted as the original administrative boundary of the city. For the
source and licence of the records the
[`crimes`](https://henriquescadbylab.github.io/spHawkes/reference/crimes.md)
events were simulated from, see
[`crimes`](https://henriquescadbylab.github.io/spHawkes/reference/crimes.md).

## See also

[`crimes`](https://henriquescadbylab.github.io/spHawkes/reference/crimes.md)

## Examples

``` r
data("crimes")
data("crimes_region")

plot(sf::st_geometry(crimes_region), col = "grey95")
plot(
  sf::st_geometry(crimes),
  add = TRUE,
  pch = 16,
  cex = 0.3
)
```

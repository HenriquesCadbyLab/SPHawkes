# Combine temporal and spatial background components into the background intensity shape (without the scalar `mu` parameter)

Multiplies the trend, weekly and monthly temporal base values by the
spatial background surface, evaluated at the observed events and
summarized over the study grid.

## Usage

``` r
background_intensity(events, trend, weekly, monthly, spatial)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)).

- trend, weekly, monthly:

  `"semihawkes_temporal_base"` objects from
  [`temporal_base_trend()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_trend.md),
  [`temporal_base_weekly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_weekly.md)
  and
  [`temporal_base_monthly()`](https://henriquescadbylab.github.io/spHawkes/reference/temporal_base_monthly.md)
  respectively.

- spatial:

  A `"semihawkes_spatial_base"` object from
  [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md).

## Value

A list with elements:

- mub_events:

  Spatial background value at each event's location.

- bg_at_events:

  `trend * weekly * monthly * mub_events` at each event – the background
  intensity shape at the observed events, i.e. \\\mu(s_i,t_i)/\mu_0\\.

- mub_grid:

  Spatial background value on the full spatial grid (same shape as
  `spatial$basevalue`).

- bg_at_all:

  Single number: the space-time integral of the background intensity
  shape over the study region and period, \\\iint_S \int_0^T
  \mu(s,t)/\mu_0 \\ dt \\ ds\\, used to normalize the log-likelihood.

- marks_events:

  Calendar marks (weekday and month) for each event, as returned by
  `.calendar_marks()`; consumed downstream by
  [`em_step()`](https://henriquescadbylab.github.io/spHawkes/reference/em_step.md).

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
bgint <- background_intensity(
  ev,
  trend   = temporal_base_trend(ev),
  weekly  = temporal_base_weekly(ev),
  monthly = temporal_base_monthly(ev),
  spatial = background_spatial(ev, cell_size = 1000)
)
summary(bgint$bg_at_events)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1896  1.4648  3.1317  3.2271  4.3051 10.1196 
bgint$bg_at_all
#> [1] 23019006210
```

# Evaluate self-exciting (triggering) contributions at events and over the study domain

For every event \\i\\, adds its triggering contribution \\g_T(t_j - t_i)
\\ g_S(x_j - x_i, y_j - y_i)\\ to every other event \\j\\ (this is `0`
whenever \\t_j \le t_i\\, because `temporal_trigger` and
`spatial_trigger` vanish for non-positive/out-of-range offsets), and
separately integrates its contribution over the full spatial grid and
`[0, TT]`. The kernels are zero beyond the ends of their own grids, so
that integral only evaluates them at the time marks and background-grid
cells within their reach of the event; the rest would add exactly zero.

## Usage

``` r
trigger_evaluate(
  events,
  background,
  temporal_trigger,
  spatial_trigger,
  exclude_self = TRUE,
  nprocs = 1
)
```

## Arguments

- events:

  A `semihawkes_events` object (from
  [`prepare_events()`](https://henriquescadbylab.github.io/spHawkes/reference/prepare_events.md)).

- background:

  A `semihawkes_spatial_base` object (from
  [`background_spatial()`](https://henriquescadbylab.github.io/spHawkes/reference/background_spatial.md))
  – only its grid and `marks` are used here, to define the spatial
  domain integrated over.

- temporal_trigger:

  A `semihawkes_trigger_temporal` object (from
  [`trigger_temporal_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_temporal_init.md),
  or a later Loop-`k` re-estimate).

- spatial_trigger:

  A `semihawkes_trigger_spatial` object (from
  [`trigger_spatial_init()`](https://henriquescadbylab.github.io/spHawkes/reference/trigger_spatial_init.md),
  or a later Loop-`k` re-estimate).

- exclude_self:

  Whether to zero each event's contribution to itself; default `TRUE`.

- nprocs:

  Parallel workers: a whole number (the default, `1`, runs serially),
  `"auto"` to choose from the machine, or a cluster from
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
  which is used and left running. `"auto"` uses one worker per physical
  core but one (two on a two-core machine), at most 2 while
  `R CMD check` runs, and runs serially for fewer than 1,000 events,
  where starting the workers costs more than they save. A number larger
  than the machine's cores is reduced to the number of cores, with a
  warning. The workers are socket clusters, which work on Windows, macOS
  and Linux, and the result is the same, to the last bit, whatever their
  number. Starting them takes up to a second, so for repeated calls pass
  a cluster;
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md)
  does this itself.

## Value

A list with `at_events` (numeric vector, one value per event: the total
triggering contribution received from all other events) and `at_all`
(single number: the same contribution integrated over the whole study
region and `[0, TT]`, summed across all events – i.e.
`mytriggers.at.events.no.A` / `mytriggers.at.all.no.A` in the original
script, without the `A` scaling parameter).

## Examples

``` r
ev <- prepare_events(sim_events, sim_region, date_col = "date",
                     end_date = "2022-01-01")
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
ev <- add_adaptive_bandwidth(ev, k = 5, min_dist = 150, min_bandwidth = 300)
ev <- add_bg_integral(ev, n_mc = 500, seed = 1)
bg <- background_spatial(ev, cell_size = 1000)
trig <- trigger_evaluate(ev, bg, trigger_temporal_init(),
                         trigger_spatial_init(x_range = 2000, x_step = 250))
summary(trig$at_events)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 2.035e-09 4.704e-09 6.385e-09 9.836e-09 2.511e-08 
trig$at_all
#> [1] 67.74254
```

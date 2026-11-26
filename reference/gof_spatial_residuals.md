# Spatial raw residuals for a fitted model

Computes observed minus expected event counts on the spatial background
grid used by the fitted semi-parametric Hawkes model.

## Usage

``` r
gof_spatial_residuals(fit, start_date = NULL, end_date = NULL)
```

## Arguments

- fit:

  A fitted object of class `semihawkes_fit` returned by
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).

- start_date:

  Optional first day of the evaluation window. If `NULL`, the start of
  the fitted observation window is used.

- end_date:

  Optional last day of the evaluation window, included in it: events on
  that day are counted, and the expected count runs to the end of the
  day. If `NULL`, the last day of the fitted observation window is used.

## Value

An object of class `semihawkes_gof_spatial`: a list containing `cells`,
`region`, `start_date` and `end_date` (the first and last days of the
window), `start_day` and `end_day` (the window as the interval
`[start_day, end_day)` in days since the start of the study) and
`window_length` (`end_day - start_day`), and `cell_size`, the side of
the background grid's cells. The `cells` data frame contains `x`, `y`,
`expected`, `observed`, `residual`, `relative_residual`,
`standardized_residual`.

## Details

For each grid cell \\B\\, the raw residual is \$\$R(B) = N(B) -
\mathbb{E}\[N(B)\],\$\$ where \\N(B)\\ is the observed number of events
in the cell and \\\mathbb{E}\[N(B)\]\\ is the model-implied expected
count.

The function also returns the relative residual,
\$\$R\_{\mathrm{rel}}(B) = \frac{N(B) - \mathbb{E}\[N(B)\]}
{\mathbb{E}\[N(B)\]},\$\$ and the Poisson-standardized residual,
\$\$R\_{\mathrm{std}}(B) = \frac{N(B) - \mathbb{E}\[N(B)\]}
{\sqrt{\mathbb{E}\[N(B)\]}}.\$\$

## Examples

``` r
# \donttest{
# A quick fit to one year of the bundled simulated data: a coarse grid, a
# coarse triggering-kernel grid and two EM iterations keep it to a few
# seconds. Real analyses need finer settings and more iterations.
cache_dir <- tempfile("spHawkes_marks_")
fit <- fit_semihawkes(
  sim_events, sim_region,
  date_col  = "date",
  end_date  = "2022-01-01",
  cell_size = 1000,
  cache_dir = cache_dir,
  trigger_spatial_args = list(x_range = 2000, x_step = 250),
  trigger_window_dist  = 2000,
  n_iter  = 2,
  verbose = FALSE,
  seed    = 1
)
#> The modelled study period ends on 2021-12-29, the day after the last event. `end_date` (2022-01-01) is used only to select events and to label outputs.
res <- gof_spatial_residuals(fit)
plot(res)

# }
```

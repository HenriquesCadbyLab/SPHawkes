# Initial (parametric) spatial triggering-function guess

Builds the Loop-0 isotropic spatial triggering kernel \\g(x,y) = 1 /
(x^2/\text{scale}^2 + y^2/\text{scale}^2 + \epsilon)\\ on a regular
grid, normalized to integrate to `1` via
[`simpson.2D()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.2D.md).

## Usage

``` r
trigger_spatial_init(
  x_range = 1000,
  y_range = x_range,
  x_step = 200,
  y_step = x_step,
  scale = 40000,
  epsilon = 0.01
)
```

## Arguments

- x_range, y_range:

  Half-width of the grid along `x`/`y` (grid spans `[-range, range]`),
  in the linear units of the working CRS (metres for a UTM projection);
  default `1000` for both, i.e. a +/-1 km window. This is how far from
  an event the triggering kernel can reach, so choose it with care: see
  *Choosing the triggering reach* in
  [`fit_semihawkes()`](https://henriquescadbylab.github.io/spHawkes/reference/fit_semihawkes.md).

- x_step, y_step:

  Grid step along `x`/`y`, in the linear units of the working CRS
  (metres for a UTM projection); default `200` for both, i.e. 200 m.
  `x_range / x_step` and `y_range / y_step` must be whole numbers (5 for
  the defaults), so that each grid is centred on 0 and has the odd
  number of points
  [`simpson.2D()`](https://henriquescadbylab.github.io/spHawkes/reference/simpson.2D.md)
  needs.

- scale, epsilon:

  Kernel shape parameters: \\g(x,y) = 1 / (x^2/\text{scale}^2 +
  y^2/\text{scale}^2 + \epsilon)\\. Defaults (`40000`, `0.01`) match the
  original script.

## Value

A list of class `"semihawkes_trigger_spatial"` with `x`, `y` (grid
coordinate sequences), `basevalue` (`length(x) x length(y)` matrix,
normalized), and `fun` (a function `fun(x, y)` that bilinearly
interpolates via
[`fields::interp.surface()`](https://rdrr.io/pkg/fields/man/interp.surface.html),
returning `0` – not `NA` – outside the grid, matching the original
`excite.spatial.fun`).

## Examples

``` r
h <- trigger_spatial_init(x_range = 2000, x_step = 100)
h
#> <semihawkes_trigger_spatial>
#>   41 x 41 grid, x in [-2000, 2000], y in [-2000, 2000]
image(h$x, h$y, h$basevalue, asp = 1,
      xlab = "x offset (m)", ylab = "y offset (m)")
```

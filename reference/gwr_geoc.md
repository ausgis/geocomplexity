# geographical complexity-geographically weighted regression

geographical complexity-geographically weighted regression

## Usage

``` r
gwr_geoc(
  formula,
  data,
  gcs = NULL,
  alpha = seq(0.05, 1, 0.05),
  bw = "RMSE",
  adaptive = TRUE,
  kernel = "gaussian"
)
```

## Arguments

- formula:

  A formula of `GCGWR` model.

- data:

  An `sf` object or spatial vector object that can be converted to `sf`
  by
  [`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html).

- gcs:

  (optional) The geocomplexity matrix corresponding to each variable,
  which is calculated by default using
  [`geocd_vector()`](https://ausgis.github.io/geocomplexity/reference/geocd_vector.md).

- alpha:

  (optional) Balancing the weights of attribute similarity matrix and
  geographic distance matrix.

- bw:

  (optional) The bandwidth used in selecting models. The optimal
  bandwidth can be selected using one of three methods: `RMSE`, `AIC`,
  and `AICc`. Default will use `RMSE`.

- adaptive:

  (optional) Whether the bandwidth value is adaptive or not. Default is
  `TRUE`.

- kernel:

  (optional) Kernel function. Default is `gaussian`.

## Value

A list with GCGWR results.

- `SDF`:

  an sf tibble with coefficients, standard errors and t values

- `diagnostic`:

  goodness of fit indicators

- `args`:

  some key parameters

## Examples

``` r
# \donttest{
## The following code takes a long time to run:
econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
g = gwr_geoc(formula = Gini ~ ., data = econineq,
             alpha = 0.5, bw = "AIC", adaptive = TRUE)
g
#> Geographical Complexity-Geographically Weighted Regression Model
#> ================================================================
#>      Kernel:  gaussian
#>   Bandwidth:  16 (Nearest Neighbours) (Optimized according to AIC)
#>       Alpha:  0.5
#> 
#> Summary of Coefficient Estimates
#> --------------------------------
#> Coefficient      Min.   1st Qu.    Median   3rd Qu.      Max.
#> Intercept     -0.073     0.319     0.333     0.345     0.424
#> Induscale     -0.403    -0.281    -0.249    -0.208    -0.063
#> IT            -0.003    -0.002    -0.002    -0.002     0.001
#> Income         0.000     0.000     0.000     0.000     0.000
#> Sexrat         0.012     0.048     0.058     0.068     0.151
#> Houseown       0.002     0.003     0.003     0.003     0.003
#> Indemp        -0.000    -0.000    -0.000    -0.000    -0.000
#> Indcom         0.000     0.000     0.000     0.000     0.000
#> Hiedu          0.000     0.001     0.001     0.002     0.002
#> 
#> Diagnostic Information
#> ----------------------
#>   RSS: 0.222
#>   ENP: 47.205
#>   EDF: 285.795
#>    R2: 0.704
#> R2adj: 0.697
#>   AIC: -1458.466
#>  AICc: -1413.992
#>  RMSE: 0.026
#> 
# }
```

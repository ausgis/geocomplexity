# constructing spatial weight matrix based on geocomplexity with spatial dependence

constructing spatial weight matrix based on geocomplexity with spatial
dependence

## Usage

``` r
geocd_swm(sfj, wt = NULL, style = "B", ...)
```

## Arguments

- sfj:

  An `sf` object or spatial vector object that can be converted to `sf`
  by
  [`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html).

- wt:

  (optional) Spatial weight matrix based on spatial adjacency or spatial
  distance relationships.

- style:

  (optional) A character that can be `B`,`W`,`C`. More to see
  [`spdep::nb2mat()`](https://r-spatial.github.io/spdep/reference/nb2mat.html).
  Default is `B`.

- ...:

  (optional) Other parameters passed to
  [`geocomplexity::geocd_vector()`](https://ausgis.github.io/geocomplexity/reference/geocd_vector.md).

## Value

A matrix

## Examples

``` r
econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
wt_gc = geocd_swm(econineq)
wt_gc[1:5,1:5]
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.0000000 0.5050696 0.4544196 0.5640641 0.5800917
#> [2,] 0.4990206 0.0000000 0.8393764 0.0000000 0.0000000
#> [3,] 0.4314267 0.8366975 0.0000000 0.0000000 0.0000000
#> [4,] 0.5440211 0.0000000 0.0000000 0.0000000 0.9429423
#> [5,] 0.5713514 0.0000000 0.0000000 0.9402259 0.0000000
```

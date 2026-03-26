# constructing spatial weight matrix based on geocomplexity with similar geographical configurations

constructing spatial weight matrix based on geocomplexity with similar
geographical configurations

## Usage

``` r
geocs_swm(sfj, wt = NULL, style = "B", ...)
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
  [`geocomplexity::geocs_vector()`](https://ausgis.github.io/geocomplexity/reference/geocs_vector.md).

## Value

A matrix

## Examples

``` r
econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
wt_gc = geocs_swm(econineq)
wt_gc[1:5,1:5]
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.6926861 0.7649766 0.7744006 0.6888728 0.7255376
#> [2,] 0.6926861 0.7649766 0.7744006 0.6888728 0.7255376
#> [3,] 0.6926861 0.7649766 0.7744006 0.6888728 0.7255376
#> [4,] 0.6926861 0.7649766 0.7744006 0.6888728 0.7255376
#> [5,] 0.6926861 0.7649766 0.7744006 0.6888728 0.7255376
```

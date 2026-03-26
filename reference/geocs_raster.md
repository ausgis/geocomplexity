# geocomplexity for spatial raster data based on geographical similarity

This function calculates geocomplexity for spatial raster data based on
geographical similarity.

## Usage

``` r
geocs_raster(r, order = 1, normalize = TRUE, similarity = 1, method = "spvar")
```

## Arguments

- r:

  `SpatRaster` object or can be converted to `SpatRaster` by
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).

- order:

  (optional) The order of the adjacency object. Default is `1`.

- normalize:

  (optional) Whether to further normalizes the calculated geocomplexity.
  Default is `TRUE`.

- similarity:

  (optional) When `similarity` is `1`, the similarity is calculated
  using geographical configuration similarity, otherwise the cosine
  similarity is calculated. Default is `1`.

- method:

  (optional) When `method` is `spvar`, variation of the similarity
  vector is represented using spatial variance, otherwise shannon
  information entropy is used. Default is `spvar`.

## Value

A `SpatRaster` object

## Note

In contrast to the
[`geocs_vector()`](https://ausgis.github.io/geocomplexity/reference/geocs_vector.md)
function, the `geocs_raster()` performs operations internally on raster
data without providing additional wt object.

## Examples

``` r
library(terra)
m1 = matrix(c(3,3,3,3,1,3,
              3,3,3,2,1,2,
              3,3,3,1,2,1,
              1,3,2,2,2,2,
              2,2,2,1,1,2,
              1,2,1,1,1,1),
           nrow = 6,
           byrow = TRUE)
m1 = rast(m1)
names(m1) = 'sim1'
m2 = m1
set.seed(123456789)
values(m2) = values(m1) + runif(ncell(m1),-1,1)
names(m2) = 'sim2'
m = c(m1,m2)
gc1 = geocs_raster(m,1)
gc2 = geocs_raster(m,2)
gc1
#> class       : SpatRaster 
#> size        : 6, 6, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 6, 0, 6  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> name        : GC 
#> min value   :  0 
#> max value   :  1 
plot(gc1)

gc2
#> class       : SpatRaster 
#> size        : 6, 6, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 6, 0, 6  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> name        : GC 
#> min value   :  0 
#> max value   :  1 
plot(gc2)

```

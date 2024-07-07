#' @title calculates geocomplexity for raster data
#' @description
#' This function calculates geocomplexity, a geospatial local complexity indicator,
#' for variable in raster data. The resulting output is A SpatRaster object.
#'
#' @param r Raster data that can be converted to `SpatRaster` by `terra::rast()`.
#'
#' @return A SpatRaster object
#' @export
#'
#' @examples
#'
#' library(terra)
#'
#' m = matrix(
#'   c(3,3,3,3,1,3,
#'     3,3,3,2,1,2,
#'     3,3,3,1,2,1,
#'     1,3,2,2,2,2,
#'     2,2,2,1,1,2,
#'     1,2,1,1,1,1),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#'
#' m = rast(m)
#' plot(m, col = c("#d2eaac", "#a3dae1", "#8cc1e1"))
#' a = geoc_raster(m)
#' a
#'
geoc_raster = \(r){
  if (!inherits(r,'SpatRaster')){
    r = terra::rast(r)
  }
  seq(1,terra::nlyr(r)) %>%
    purrr::map(\(i) terra::app(r[[i]],standardize_vector)) %>%
    terra::rast() -> r
  geocres = terra::focalCpp(r,w = 3,
                            RasterGeocSpatialDependence,
                            fillvalue = NA)
  seq(1,terra::nlyr(geocres)) %>%
    purrr::map(\(i) terra::app(geocres[[i]],normalize_vector)) %>%
    terra::rast() -> geocres
  return(geocres)
}

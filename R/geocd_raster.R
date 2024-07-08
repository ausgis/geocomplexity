#' @title calculates geocomplexity for raster data based on spatial dependence
#' @description
#' This function calculates geocomplexity for raster data based on spatial dependence.
#'
#' @param r Raster object that can be converted to `SpatRaster` by `terra::rast()`.
#' @param order (optional) The order of the adjacency object. Default is `1`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#'
#' @return A SpatRaster object
#' @export
#'
#' @examples
#' library(terra)
#' m = matrix(c(3,3,3,3,1,3,
#'              3,3,3,2,1,2,
#'              3,3,3,1,2,1,
#'              1,3,2,2,2,2,
#'              2,2,2,1,1,2,
#'              1,2,1,1,1,1),
#'            nrow = 6,
#'            byrow = TRUE)
#' m = rast(m)
#' names(m) = 'sim'
#' plot(m, col = c("#d2eaac", "#a3dae1", "#8cc1e1"))
#' gc1 = geocd_raster(m,1)
#' gc2 = geocd_raster(m,2)
#' gc1
#' gc2
#'
geocd_raster = \(r,order = 1,normalize = TRUE){
  if (!inherits(r,'SpatRaster')){
    r = terra::rast(r)
  }
  rastlayername = names(r)
  seq(1,terra::nlyr(r)) %>%
    purrr::map(\(i) terra::app(r[[i]],standardize_vector)) %>%
    terra::rast() -> r
  geocres = terra::focalCpp(r, w = 2*order + 1,
                            RasterGeoCDependence,
                            fillvalue = NA)
  if (normalize) {
    seq(1,terra::nlyr(geocres)) %>%
      purrr::map(\(i) terra::app(geocres[[i]],normalize_vector)) %>%
      terra::rast() -> geocres
  }
  names(geocres) = paste0('Geocomplexity_',rastlayername)
  return(geocres)
}

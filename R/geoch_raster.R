#' @title geocomplexity for raster data based on spatial stratified heterogeneity
#' @description
#' This function calculates geocomplexity for raster data based on spatial stratified heterogeneity.
#'
#' @param r Raster object that can be converted to `SpatRaster` by `terra::rast()`.
#' @param order (optional) The order of the adjacency object. Default is `1`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param method (optional) When `method` is `spvar`, variation of the attribute vector is
#' represented using spatial variance, otherwise Shannon information entropy is used. Default
#' is `spvar`.
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
#' gc1 = geoch_raster(m,1)
#' gc2 = geoch_raster(m,2)
#' gc1
#' gc2
#'
geoch_raster = \(r,order = 1,normalize = TRUE,method = 'spvar'){
  if (!inherits(r,'SpatRaster')){
    r = terra::rast(r)
  }
  rastlayername = names(r)
  seq(1,terra::nlyr(r)) %>%
    purrr::map(\(i) terra::app(r[[i]],normalize_vector)) %>%
    terra::rast() -> r
  rmat = as.matrix(r)
  imat = seq(0,terra::ncell(r[[1]])-1) %>%
    as.integer() %>%
    matrix(nrow = terra::nrow(r[[1]]), byrow = TRUE)
  geocres = seq(1,terra::nlyr(r)) %>%
    purrr::map(\(.i) r[[.i]] %>%
                 terra::values() %>%
                 RasterGeoCSSH(x = ., iw = imat,
                               w = as.integer(2*order+1),
                               method = method)) %>%
    terra::rast()
  if (normalize) {
    seq(1,terra::nlyr(geocres)) %>%
      purrr::map(\(i) terra::app(geocres[[i]],normalize_vector)) %>%
      terra::rast() -> geocres
  }
  names(geocres) = paste0('Geocomplexity_',rastlayername)
  return(geocres)
}

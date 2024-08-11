#' @title geocomplexity for raster data based on geographical similarity
#' @description
#' This function calculates geocomplexity for raster data based on geographical similarity.
#' @note
#' In contrast to the `geocs_vector()` function, the `geocs_raster()` performs operations
#' internally on raster data without providing additional wt object.
#'
#' @param r Raster object that can be converted to `SpatRaster` by `terra::rast()`.
#' @param order (optional) The order of the adjacency object. Default is `1`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param similarity (optional) When `similarity` is `1`, the similarity is calculated using
#' geographical configuration similarity, otherwise the cosine similarity is calculated.
#' Default is `1`.
#' @param method (optional) When `method` is `spvar`, variation of the similarity vector is
#' represented using spatial variance, otherwise shannon information entropy is used. Default
#' is `spvar`.
#'
#' @return A SpatRaster object
#' @export
#'
#' @examples
#' library(terra)
#' m1 = matrix(c(3,3,3,3,1,3,
#'               3,3,3,2,1,2,
#'               3,3,3,1,2,1,
#'               1,3,2,2,2,2,
#'               2,2,2,1,1,2,
#'               1,2,1,1,1,1),
#'            nrow = 6,
#'            byrow = TRUE)
#' m1 = rast(m1)
#' names(m1) = 'sim1'
#' m2 = m1
#' set.seed(123456789)
#' values(m2) = values(m1) + runif(ncell(m1),-1,1)
#' names(m2) = 'sim2'
#' m = c(m1,m2)
#' gc1 = geocs_raster(m,1)
#' gc2 = geocs_raster(m,2)
#' gc1
#' gc2
#'
geocs_raster = \(r, order = 1, normalize = TRUE,
                 similarity = 1,method = 'spvar'){
  if (!inherits(r,'SpatRaster')){
    r = terra::rast(r)
  }
  rastlayername = names(r) %>%
    paste0(collapse = '_')
  seq(1,terra::nlyr(r)) %>%
    purrr::map(\(i) terra::app(r[[i]],sdsfun::normalize_vector)) %>%
    terra::rast() -> r
  if (terra::nlyr(r) == 1) {
    stop('To use `geocs_raster`, the number of layers in r must be greater than or equal to 2')
  }
  rmat = as.matrix(r)
  imat = seq(0,terra::ncell(r[[1]])-1) %>%
    as.integer() %>%
    matrix(nrow = terra::nrow(r[[1]]), byrow = TRUE)
  geocres = r[[1]]
  geocres = RasterGeoCSimilarity(rmat,imat,as.integer(2*order+1),similarity,method)
  if (normalize) {
    geocres = sdsfun::normalize_vector(geocres)
  }
  r1 = r[[1]]
  terra::values(r1) = geocres
  names(r1) = paste0('Geocomplexity_',rastlayername)
  return(r1)
}

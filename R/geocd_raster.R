#' @title geocomplexity for spatial raster data based on spatial dependence
#' @description
#' This function calculates geocomplexity for spatial raster data based on spatial dependence.
#' @references
#' Zehua Zhang, Yongze Song, Peng Luo & Peng Wu (2023) Geocomplexity explains spatial errors,
#' International Journal of Geographical Information Science, 37:7, 1449-1469,
#' DOI: 10.1080/13658816.2023.2203212
#'
#' Anselin, L. (2019). A local indicator of multivariate spatial association: Extending
#' geary’s c. Geographical Analysis, 51(2), 133–150. https://doi.org/10.1111/gean.12164
#'
#' @note
#' In contrast to the `geocd_vector()` function, the `geocd_raster()` performs operations
#' internally on raster data based on neighborhood operations(focal) without providing
#' additional wt object.
#'
#' @param r `SpatRaster` object or can be converted to `SpatRaster` by `terra::rast()`.
#' @param order (optional) The order of the adjacency object. Default is `1`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param method (optional) In instances where the method is `moran`, geocomplexity is
#' determined using local moran measure method. Conversely, when the method is `spvar`,
#' the spatial variance of local attribute data serves to characterize geocomplexity. For all
#' other methods, the shannon information entropy of attribute data is employed to represent
#' geocomplexity. The selection of the method can be made from any one of the three options:
#' `moran`, `spvar` or `entropy`. Default is `moran`.
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
#' plot(gc1)
#' gc2
#' plot(gc2)
#'
geocd_raster = \(r,order = 1,normalize = TRUE,method = 'moran'){
  if (!inherits(r,'SpatRaster')){
    r = terra::rast(r)
  }
  rastlayername = names(r)
  seq(1,terra::nlyr(r)) %>%
    purrr::map(\(i) terra::app(r[[i]],sdsfun::standardize_vector)) %>%
    terra::rast() -> r

  if (method == 'moran'){
    geocres = terra::focalCpp(r, w = 2*order + 1,
                              RasterGeoCMoran,
                              fillvalue = NA)
  } else {
    imat = seq(0,terra::ncell(r[[1]])-1) %>%
      as.integer() %>%
      matrix(nrow = terra::nrow(r[[1]]), byrow = TRUE)
    geocres = r
    for (.i in seq(1,terra::nlyr(r))) {
      terra::values(geocres[[.i]]) = r[[.i]] %>%
        terra::values() %>%
        RasterGeoCSSH(x = ., iw = imat,
                      w = as.integer(2*order+1),
                      method = method)
    }
  }

  if (normalize) {
    seq(1,terra::nlyr(geocres)) %>%
      purrr::map(\(i) terra::app(geocres[[i]],sdsfun::normalize_vector)) %>%
      terra::rast() -> geocres
  }

  names(geocres) = paste0("GC_",rastlayername)
  return(geocres)
}

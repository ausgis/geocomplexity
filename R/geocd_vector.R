#' @title geocomplexity for vector data based on spatial dependence
#' @description
#' This function calculates geocomplexity for vector data based on spatial dependence.
#' @references
#' Zehua Zhang, Yongze Song, Peng Luo & Peng Wu (2023) Geocomplexity explains spatial errors,
#' International Journal of Geographical Information Science, 37:7, 1449-1469,
#' DOI: 10.1080/13658816.2023.2203212
#'
#' Anselin, L. (2019). A local indicator of multivariate spatial association: Extending
#' geary’s c. Geographical Analysis, 51(2), 133–150. https://doi.org/10.1111/gean.12164
#'
#' @details
#' The formula for geocomplexity which uses local moran measure method is
#'
#' \eqn{\rho_i = -\frac{1}{m} Z_i \sum\limits_{j=1}^m W_{ij} Z_j -\frac{1}{m} \sum\limits_{j=1}^m W_{ij} Z_j \frac{1}{V_{k}}\sum\limits_{k=1}^n W_{jk} W_{ik} Z_k}
#'
#' @note
#' If `wt` is not provided, for polygon vector data, `geocomplexity` will use a first-order queen
#' adjacency binary matrix; for point vector data, the six nearest points are used as adjacency
#' objects to generate an adjacency binary matrix.
#'
#' @param sfj An `sf` object or vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param method (optional) In instances where the method is `moran`, geocomplexity is
#' determined using local moran measure method. Conversely, when the method is `spvar`,
#' the spatial variance of attribute data serves to characterize geocomplexity. For all
#' other methods, the shannon information entropy of attribute data is employed to represent
#' geocomplexity. The selection of the method can be made from any one of the three options:
#' `moran`, `spvar` or `entropy`. Default is `moran`.
#'
#' @return An sf object
#' @export
#'
#' @examples
#' data("income")
#' inc = dplyr::select(income,Income)
#' gc = geocd_vector(inc)
#' gc
#'
#' library(ggplot2)
#' library(viridis)
#' ggplot(gc) +
#'    geom_sf(aes(fill = Geocomplexity_Income)) +
#'    scale_fill_viridis(option="mako", direction = -1) +
#'    theme_bw()
#'
geocd_vector = \(sfj,wt = NULL,normalize = TRUE,method = 'moran'){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  if (is.null(wt)){
    if (sdsfun::sf_geometry_type(sfj) %in% c('polygon','multipolygon')){
      wt = sdsfun::spdep_contiguity_swm(sfj,
                                        queen = TRUE,
                                        style = 'B',
                                        zero.policy = TRUE)
    } else if (sdsfun::sf_geometry_type(sfj) %in% c('point','multipoint')){
      wt = sdsfun::spdep_contiguity_swm(sfj,
                                        k = 6,
                                        style = 'B',
                                        zero.policy = TRUE)
    } else {
      stop('Only support (multi-)point or (multi-)polygon vector data!')
    }
  } else {
    wt = check_wt(wt)
  }
  sfj_attr = sf::st_drop_geometry(sfj) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                sdsfun::standardize_vector))
  vectlayername = names(sfj_attr)

  if (method == 'moran'){
    geocvec = dplyr::mutate(sfj_attr,
                            dplyr::across(dplyr::everything(),
                                          \(.x) VectorGeoCMoran(.x,wt)))
  } else {
    geocvec = dplyr::mutate(sfj_attr,
                            dplyr::across(dplyr::everything(),
                                          \(.x) VectorGeoCSSH(.x,wt,method)))
  }

  if (normalize) {
    geocvec = dplyr::mutate(geocvec,
                            dplyr::across(dplyr::everything(),
                                          sdsfun::normalize_vector))
  }

  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

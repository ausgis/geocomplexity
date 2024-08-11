#' @title geocomplexity for vector data based on geographical similarity
#' @description
#' This function calculates geocomplexity for in vector data based on geographical similarity.
#' @details
#' The geographical configuration similarity is calculated as:
#'
#' \eqn{S({\bf u}_{\alpha},{\bf v}_{\beta})=min\{E_{i}(e_{i}({\bf u}_{\alpha}),e_{i}({\bf v}_{\beta}))\}}
#'
#' \eqn{E_{i}({\bf u}_{\alpha},{\bf v}_{\beta})=\exp\left(-{\frac{\left(e_{i}({\bf u}_{\alpha})-e_{i}({\bf v}_{\beta})\right)^{2}}{2\left(\sigma^{2}/\delta({\bf v}_{\beta})\right)^{2}}}\right)}
#'
#' \eqn{\delta({\bf u}_{\alpha},{\bf v})=\sqrt{\frac{\sum_{\beta=1}^{n}(e({\bf u}_{\alpha})-e({\bf v}_{\beta}))^{2}}{n}}}
#'
#' The spatial variance is calculated as:
#'
#' \eqn{\Gamma = \frac{\sum_i \sum_{j \neq i} \omega_{ij}\frac{(y_i-y_j)^2}{2}}{\sum_i \sum_{j \neq i} \omega_{ij}}}
#'
#' @param sfj An `sf` object or vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. If `wt` is not
#' provided, `geocomplexity` will use a first-order inverse distance weight matrix via
#' `sdsfun::inverse_distance_swm()` function.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param similarity (optional) When `similarity` is `1`, the similarity is calculated using
#' geographical configuration similarity, otherwise the cosine similarity is calculated.
#' Default is `1`.
#' @param method (optional) When `method` is `spvar`, variation of the similarity vector is
#' represented using spatial variance, otherwise shannon information entropy is used. Default
#' is `spvar`.
#'
#' @return An sf object
#' @export
#'
#' @examples
#' data("income")
#' income = dplyr::select(income,-SA3_CODE16)
#' gc = geocs_vector(income)
#' gc
#'
#' library(ggplot2)
#' library(viridis)
#' ggplot(gc) +
#'    geom_sf(aes(fill = Geocomplexity_Income_Gini)) +
#'    scale_fill_viridis(option="mako", direction = -1) +
#'    theme_bw()
#'
geocs_vector = \(sfj, wt = NULL, normalize = TRUE,
                 similarity = 1, method = 'spvar'){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  if (is.null(wt)){
    wt = sdsfun::inverse_distance_swm(sfj)
  } else {
    wt = check_wt(wt)
  }
  sfj_attr = sf::st_drop_geometry(sfj) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                sdsfun::normalize_vector))
  if (ncol(sfj_attr) == 1) {
    stop('To use `geocs_vector`, the number of attribute columns in sfj must be greater than or equal to 2')
  }
  vectlayername = paste(names(sfj_attr),collapse = '_')
  geocvec = VectorGeoCSimilarity(as.matrix(sfj_attr),wt,similarity,method)
  if (normalize) {
    geocvec = sdsfun::normalize_vector(geocvec)
  }
  geocvec = tibble::as_tibble_col(geocvec)
  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

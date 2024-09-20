#' @title geocomplexity for spatial vector data based on geographical similarity
#' @description
#' This function calculates geocomplexity for in spatial vector data based on geographical similarity.
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
#' @param sfj An `sf` object or spatial vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. If `wt` is not
#' provided, `geocomplexity` will use a first-order inverse distance weight matrix via
#' `sdsfun::inverse_distance_swm()` function.
#' @param method (optional) When `method` is `spvar`, variation of the similarity vector is
#' represented using spatial variance, otherwise shannon information entropy is used. Default
#' is `spvar`.
#' @param similarity (optional) When `similarity` is `1`, the similarity is calculated using
#' geographical configuration similarity, otherwise the cosine similarity is calculated.
#' Default is `1`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param returnsf (optional) When `returnsf` is `TRUE`, return an sf object, otherwise a tibble.
#' Default is `TRUE`.
#'
#' @return A tibble (`returnsf` is `FALSE`) or an sf object (`returnsf` is `TRUE`)
#' @export
#'
#' @examples
#' econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
#' gc = geocs_vector(dplyr::select(econineq,-Gini))
#' gc
#'
#' library(ggplot2)
#' library(viridis)
#' ggplot(gc) +
#'    geom_sf(aes(fill = GC)) +
#'    scale_fill_viridis(option="mako", direction = -1) +
#'    theme_bw()
#'
geocs_vector = \(sfj, wt = NULL,
                 method = 'spvar', similarity = 1,
                 normalize = TRUE, returnsf = TRUE){
  if (!inherits(sfj,'sf')) {
    sfj = sf::st_as_sf(sfj)
  }

  if (is.null(wt)) {
    wt = sdsfun::inverse_distance_swm(sfj)
  }

  sfj_attr = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                sdsfun::normalize_vector))
  if (ncol(sfj_attr) == 1) {
    stop('To use `geocs_vector`, the number of attribute columns in sfj must be greater than or equal to 2')
  }
  geocvec = VectorGeoCSimilarity(as.matrix(sfj_attr),wt,similarity,method)

  if (normalize) {
    geocvec = sdsfun::normalize_vector(geocvec)
  }

  geocvec = tibble::as_tibble_col(geocvec)
  names(geocvec) = "GC"

  if (returnsf) {
    geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  }

  return(geocvec)
}

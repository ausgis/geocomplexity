#' @title geocomplexity for vector data based on geographical similarity
#' @description
#' This function calculates geocomplexity for in vector data based on geographical similarity.
#'
#' @param sfj Vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order inverse distance weight matrix via
#' `inverse_distance_weight()`.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#'
#' @return An sf object
#' @export
#'
#' @examples
#' library(sf)
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
geocs_vector = \(sfj,wt = NULL,normalize = TRUE){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  if (is.null(wt)){
    wt = inverse_distance_weight(sfj)
  }
  sfj_attr = sf::st_drop_geometry(sfj) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                standardize_vector))
  if (ncol(sfj_attr) == 1) {
    stop('To use `geocs_vector`, the number of attribute columns in sfj must be greater than or equal to 2')
  }
  vectlayername = paste(names(sfj_attr),collapse = '_')
  geocvec = VectorGeoCSimilarity(as.matrix(sfj_attr),wt)
  if (normalize) {
    geocvec = normalize_vector(geocvec)
  }
  geocvec = tibble::as_tibble_col(geocvec)
  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

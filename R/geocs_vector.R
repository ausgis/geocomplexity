#' @title calculates geocomplexity for vector data based on spatial similarity
#' @description
#' This function calculates geocomplexity, a geospatial local complexity indicator,
#' for variable in vector data based on spatial similarity.
#'
#' @param sfj Vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order inverse distance weight matrix via
#' `inverse_distance_weight()`.
#' @param normalize (optional) Whether to further normalizes the spatial local complexity.
#' Default is `TRUE`.
#'
#' @return An sf object
#' @export
#'
#' @examples
#' library(sf)
#' data("income")
#' inc = dplyr::select(income,Income)
#' gc = geocs_vector(inc)
#' gc
#'
#' library(ggplot2)
#' library(viridis)
#' ggplot(gc) +
#'    geom_sf(aes(fill = Geocomplexity_Income)) +
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
  vectlayername = names(sfj_attr)
  geocvec = dplyr::mutate(sfj_attr,
                          dplyr::across(dplyr::everything(),
                                        \(.x) VectorGeoCSimilarity(.x,wt)))
  if (normalize) {
    geocvec = dplyr::mutate(geocvec,
                            dplyr::across(dplyr::everything(),
                                          normalize_vector))
  }
  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

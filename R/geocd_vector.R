#' @title calculates geocomplexity for vector data based on spatial dependence
#' @description
#' This function calculates geocomplexity for vector data based on spatial dependence.
#'
#' @param sfj Vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order queen adjacency binary matrix.
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
geocd_vector = \(sfj,wt = NULL,normalize = TRUE){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  if (is.null(wt)){
    nb_queen = spdep::poly2nb(sfj, queen=TRUE)
    wt = spdep::nb2mat(nb_queen, style='B',
                       zero.policy = TRUE)
  }
  sfj_attr = sf::st_drop_geometry(sfj) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                standardize_vector))
  vectlayername = names(sfj_attr)
  geocvec = dplyr::mutate(sfj_attr,
                          dplyr::across(dplyr::everything(),
                                        \(.x) VectorGeoCDependence(.x,wt)))
  if (normalize) {
    geocvec = dplyr::mutate(geocvec,
                            dplyr::across(dplyr::everything(),
                                          normalize_vector))
  }
  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

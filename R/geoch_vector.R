#' @title geocomplexity for vector data based on spatial stratified heterogeneity
#' @description
#' This function calculates geocomplexity for in vector data based on spatial stratified heterogeneity.
#'
#' @param sfj An `sf` object or vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order queen adjacency binary matrix via
#' `spdep` package.
#' @param normalize (optional) Whether to further normalizes the calculated geocomplexity.
#' Default is `TRUE`.
#' @param method (optional) When `method` is `spvar`, variation of the attribute vector is
#' represented using spatial variance, otherwise Shannon information entropy is used. Default
#' is `spvar`.
#'
#' @return An sf object
#' @export
#'
#' @examples
#' data("income")
#' inc = dplyr::select(income,Income)
#' gc = geoch_vector(inc)
#' gc
#'
#' library(ggplot2)
#' library(viridis)
#' ggplot(gc) +
#'    geom_sf(aes(fill = Geocomplexity_Income)) +
#'    scale_fill_viridis(option="mako", direction = -1) +
#'    theme_bw()
#'
geoch_vector = \(sfj,wt = NULL,normalize = TRUE,method = 'spvar'){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  if (is.null(wt)){
    nb_queen = spdep::poly2nb(sfj, queen=TRUE)
    wt = spdep::nb2mat(nb_queen, style='B',
                       zero.policy = TRUE)
  } else {
    wt = check_wt(wt)
  }
  sfj_attr = sf::st_drop_geometry(sfj) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                normalize_vector))
  vectlayername = names(sfj_attr)
  geocvec = dplyr::mutate(sfj_attr,
                          dplyr::across(dplyr::everything(),
                                        \(.x) VectorGeoCSSH(.x,wt,method)))
  if (normalize) {
    geocvec = dplyr::mutate(geocvec,
                            dplyr::across(dplyr::everything(),
                                          normalize_vector))
  }
  names(geocvec) = paste0('Geocomplexity_',vectlayername)
  geocvec = sf::st_set_geometry(geocvec,sf::st_geometry(sfj))
  return(geocvec)
}

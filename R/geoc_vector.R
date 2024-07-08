geoc_vector = \(sfj,wt = NULL,normalize = TRUE){
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
  geocvec = dplyr::mutate(sfj_attr,
                          dplyr::across(dplyr::everything(),
                                        \(.x) VectorGeoCDependence(.x,wt)))
  if (normalize) {
    geocvec = dplyr::mutate(geocvec,
                            dplyr::across(dplyr::everything(),
                                          normalize_vector))
  }
  return(geocvec)
}

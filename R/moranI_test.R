moran_test = \(sfj,wt = NULL,
               alternative = "greater",
               symmetrize = TRUE){
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
  dmat = sfj %>%
    sf::st_drop_geometry() %>%
    as.matrix()
  if (!is.null(colnames(dmat))) {
    xnams = colnames(dmat)
  }
  if (0 %in% apply(dmat, 2, sd)) {
    warning("Constant term detected in attribute columns for `sfj`")
  }

  if (anyNA(dmat) | anyNA(wt)) {
    stop("Missing values detected")
  }
  if (!(alternative %in% c("greater", "lower", "two.sided"))) {
    stop("Invalid input: 'alternative' must be either 'greater',
         'lower', or 'two.sided'")
  }

}

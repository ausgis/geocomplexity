#' @title global spatial autocorrelation test
#' @description
#' Spatial autocorrelation test based on global Moran'I index.
#' @note
#' This is a `C++` implementation of the `MI.vec` function in `spfilteR` package,
#' and embellishes the console output.
#'
#' The return result of this function is actually a `list`, please access the result
#' tibble using `$result`.
#'
#' @param sfj An `sf` object or vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order queen adjacency binary matrix via
#' `spdep` package.
#' @param alternative
#' @param symmetrize
#'
#' @return
#' @export
#'
#' @examples
moran_test = \(sfj, wt = NULL,
               alternative = "greater",
               symmetrize = FALSE){
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
  if (!(alternative %in% c("greater", "lower", "two.sided"))) {
    stop("Invalid input: `alternative` must be either `greater`,
         `lower`, or `two.sided`")
  }
  mitres = MI_vec(dmat, wt, alternative, symmetrize)
  mitres = tibble::as_tibble(mitres) %>%
    dplyr::mutate(Variable = colnames(dmat)) %>%
    dplyr::select(dplyr::all_of("Variable"),
                  dplyr::everything())
  res = list(result = mitres)
  class(res) = 'moran_test'
  return(res)
}

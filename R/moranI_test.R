#' @title global spatial autocorrelation test
#' @description
#' Spatial autocorrelation test based on global Moran'I Index.
#' @note
#' This is a `C++` implementation of the `MI.vec` function in `spfilteR` package,
#' and embellishes the console output.
#'
#' The return result of this function is actually a `list`, please access the result
#' tibble using `$result`.
#'
#' The non-numeric columns of the attribute columns in `sfj` are ignored.
#'
#' @param sfj An `sf` object or vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix. Must be a `matrix` class. You can get a
#' spatial weight matrix from `spdep`,`rgeoda` or `tidyrgeoda` package. If `wt` is not
#' provided, `geocomplexity` will use a first-order queen adjacency binary matrix via
#' `spdep` package.
#' @param alternative (optional) Specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'.
#' @param symmetrize (optional) Whether or not to symmetrize the asymmetrical spatial weight matrix
#' \emph{\strong{wt}} by: 1/2 * (\emph{\strong{wt}} + \emph{\strong{wt}}'). Default is `FALSE`.
#'
#' @return A list with `moran_test` class and result stored on the `result` tibble.
#' Which contains the following information for each variable:
#' \describe{
#' \item{\code{I}}{observed value of the Moran coefficient}
#' \item{\code{EI}}{expected value of Moran's I}
#' \item{\code{VarI}}{variance of Moran's I (under normality)}
#' \item{\code{zI}}{standardized Moran coefficient}
#' \item{\code{pI}}{\emph{p}-value of the test statistic}
#' }
#'
#' @export
#'
#' @examples
#' data("income")
#' moran_test(income)
#'
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
  if (!(alternative %in% c("greater", "less", "two.sided"))) {
    stop("Invalid input: `alternative` must be either `greater`, `less`, or `two.sided`")
  }
  dmat = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    as.matrix()
  mitres = MI_vec(dmat, wt, alternative, symmetrize)
  mitres = tibble::as_tibble(mitres) %>%
    dplyr::mutate(Variable = colnames(dmat)) %>%
    dplyr::select(dplyr::all_of("Variable"),
                  dplyr::everything())
  res = list(result = mitres)
  class(res) = 'moran_test'
  return(res)
}

#' print global spatial autocorrelation test result
#' @export
#' @noRd
print.moran_test = \(x,...){
  # cat("\n ***               global spatial autocorrelation test               ")
  # cat("\n ---------------------------------------------------------------------")
  # pander::pander(x$result)
  x = as.data.frame(x$result)
  print_global_moranI(x)
}

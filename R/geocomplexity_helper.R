#' @title standardization
#'
#' @param x A numeric vector
#'
#' @return A standardized numeric vector
#' @export
#'
#' @examples
#' standardize_vector(1:10)
#'
standardize_vector = \(x){
  return((x - mean(x,na.rm = TRUE)) / stats::sd(x,na.rm = TRUE))
}

#' @title normalization
#'
#' @param x A continuous numeric vector.
#'
#' @return A continuous vector which has normalized.
#' @export
#'
#' @examples
#' normalize_vector(c(-5,1,5,0.01,0.99))
#'
normalize_vector = \(x){
  xmin = range(x,na.rm = TRUE)[1]
  xmax = range(x,na.rm = TRUE)[2]
  xnew = (x - xmin) / (xmax - xmin)
  return(xnew)
}

#' @title calculate inverse distance weight
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#' @description
#' Function for calculate inverse distance weight.
#' @details
#' The inverse distance weight formula is
#' \eqn{w_{ij} = 1 / d_{ij}^\alpha}
#'
#' @param locx The x axis location.
#' @param locy The y axis location.
#' @param power (optional) Default is 1. Set to 2 for gravity weights.
#' @param is_arc (optional) FALSE (default) or TRUE, whether to compute arc distance.
#'
#' @return A inverse distance weight matrices with class of `matrix`.
#' @export
#'
#' @examples
#' x = 1:10
#' y = 1:10
#' inverse_distance_weight(x,y)
#' inverse_distance_weight(x,y,is_arc = TRUE)
#'
inverse_distance_weight = \(locx,locy,power = 1,is_arc = FALSE){
  if (is_arc) {
    distij = stats::as.dist(geosphere::distm(matrix(c(locx,locy),
                                                    ncol = 2)))
  } else {
    distij = stats::dist(data.frame(locx,locy))
  }
  wij = 1 / distij ^ power
  return(as.matrix(wij))
}

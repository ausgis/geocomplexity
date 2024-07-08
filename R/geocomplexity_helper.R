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
#' @param sfj Vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param power (optional) Default is 1. Set to 2 for gravity weights.
#'
#' @return A inverse distance weight matrices with class of `matrix`.
#' @export
#'
#' @examples
#' data("income")
#' wt = inverse_distance_weight(income)
#' wt[1:5,1:5]
#'
inverse_distance_weight = \(sfj,power = 1){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  is_arc = ifelse(sf::st_is_longlat(sfj),TRUE,FALSE)
  coords = sfj %>%
    sf::st_point_on_surface() %>%
    sf::st_coordinates() %>%
    {.[,c('X','Y')]}
  if (is_arc) {
    distij = stats::as.dist(geosphere::distm(coords))
  } else {
    distij = stats::dist(as.data.frame(coords))
  }
  wij = 1 / distij ^ power
  return(as.matrix(wij))
}

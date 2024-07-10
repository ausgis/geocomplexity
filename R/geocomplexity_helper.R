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
  if (is_arc) {
    coords = suppressWarnings(sf::st_point_on_surface(sfj)) %>%
      sf::st_coordinates() %>%
      {.[,c('X','Y')]}
    distij = stats::as.dist(geosphere::distm(coords))
  } else {
    coords = suppressWarnings(sf::st_centroid(sfj)) %>%
      sf::st_coordinates() %>%
      {.[,c('X','Y')]}
    distij = stats::dist(as.data.frame(coords))
  }
  wij = 1 / distij ^ power
  return(as.matrix(wij))
}

#' @title get the geometry column name of an sf object
#'
#' @param sfj An `sf` object.
#'
#' @return A character.
#' @export
#'
#' @examples
#' data("income")
#' st_geometry_name(income)
#'
st_geometry_name = \(sfj){
  if (!inherits(sfj,'sf')){
    sfj = sf::st_as_sf(sfj)
  }
  gname = attr(sfj, "sf_column")
  return(gname)
}

#' @title check whether wt is a matrix class
#'
#' @param wt Spatial weight matrix
#'
#' @return If `wt` is a `matrix` class, return `wt` itself, otherwise an error is raised
#' and execution stops.
#'
#' @examples
#' data("income")
#' wt = inverse_distance_weight(income)
#' check_wt(wt)
#'
check_wt = \(wt){
  if (!any(class(wt) %in% c("matrix", "Matrix"))) {
    stop("wt must be of class `matrix`")
  }
  if (any(class(wt) != "matrix")) {
    wt = as.matrix(wt)
  }
  return(wt)
}

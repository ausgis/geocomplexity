#' @title Pipe operator
#' @description
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @title standardization
#' @description
#' To calculate the Z-score using variance normalization, the formula is as follows:
#'
#' \eqn{Z = \frac{(x - mean(x))}{sd(x)}}
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

#' @title check whether wt is a matrix class
#' @description
#' If `wt` is a `matrix` class, return `wt` itself, othe rwise an error is raised
#' and execution stops.
#' @keywords internal
#'
#' @param wt Spatial weight matrix
#'
#' @return A `matrix`.
#' @export
#'
#' @examples
#' data("income")
#' wt = sdsfun::inverse_distance_swm(income)
#' wt = check_wt(wt)
#' wt[1:5,1:5]
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

#' @title check sf geometry type
#' @keywords internal
#'
#' @param sfj
#'
#' @return A character
#' @export
#'
#' @examples
#' data("income")
#' check_geometry_type(income)
#'
check_geometry_type = \(sfj){
  sfj_type = sfj %>%
    sf::st_geometry_type() %>%
    as.character() %>%
    unique()

  if (length(sfj_type) != 1) {
    stop('Please keep one geometry type in an sf object!')
  }

  return(tolower(sfj_type))
}

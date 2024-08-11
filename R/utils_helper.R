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

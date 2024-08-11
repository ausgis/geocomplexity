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

#' check whether wt is a matrix
#' @noRd
.check_wt = \(wt){
  if (!any(class(wt) %in% c("matrix", "Matrix"))) {
    stop("wt must be of class `matrix`")
  }
  if (any(class(wt) != "matrix")) {
    wt = as.matrix(wt)
  }
  return(wt)
}

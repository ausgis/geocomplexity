#' @title spatial weight matrix based on geographical complexity
#'
#' @param gc Geographical complexity vector.
#' @param wt Spatial weight matrix based on spatial adjacency or spatial distance relationships.
#' @param style (optional) A character that can be `B`,`W`,`C`.  More to see `spdep::nb2mat()`.
#' Default is `W`.
#'
#' @return A matrix
#' @export
#'
#' @examples
#' data("income")
#' inc = dplyr::select(income,Income)
#' gc = inc %>%
#'   geocd_vector(returnsf = FALSE) %>%
#'   dplyr::pull(1)
#' wt = sdsfun::spdep_contiguity_swm(inc,style = 'B')
#' wt_gc = geoc_swm(gc,wt)
#' wt_gc[1:5,1:5]
#'
geoc_swm = \(gc,wt,style = 'W'){
  wt_gc = GeoCSWT(gc,wt,style)
  return(wt_gc)
}

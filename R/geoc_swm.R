#' @title spatial weight matrix based on geographical complexity
#'
#' @param sfj An `sf` object or spatial vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param wt (optional) Spatial weight matrix based on spatial adjacency or spatial distance relationships.
#' @param style (optional) A character that can be `B`,`W`,`C`.  More to see `spdep::nb2mat()`.
#' Default is `B`.
#' @param ... (optional) Other parameters passed to `geocomplexity::geocs_vector()`.
#'
#' @return A matrix
#' @export
#'
#' @examples
#' econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
#' wt_gc = geocs_swm(econineq)
#' wt_gc[1:5,1:5]
#'
geocs_swm = \(sfj, wt = NULL, style = 'B', ...){
  if (is.null(wt)) {
    wt = sdsfun::inverse_distance_swm(sfj)
  }

  gc = dplyr::pull(geocomplexity::geocs_vector(sfj,wt,...),1)
  wt_gc = GeoCS_SWM(gc,wt,style)

  return(wt_gc)
}

#' Geographical Complexity-Geographically Weighted Regression
#'
#' @param formula A formula of `GCGWR` model.
#' @param data An `sf` object or spatial vector object that can be converted to `sf` by `sf::st_as_sf()`.
#' @param gcs (optional) The geocomplexity matrix corresponding to each variable, which is calculated
#' by default using `geocd_vector()`.
#' @param alpha (optional) Balancing the weights of attribute similarity matrix and geographic distance matrix.
#' @param bw (optional) The bandwidth used in selecting models. The optimal bandwidth will be
#' selected based on `RMSE` by default.
#' @param adaptive (optional) Whether the bandwidth value is adaptive or not. Default is `TRUE`.
#' @param kernel (optional) Kernel function. Default is `gaussian`.
#'
#' @return A list with GCGWR results.
#' \describe{
#' \item{\code{SDF}}{an sf tibble with coefficients, standard errors and t values}
#' \item{\code{diagnostic}}{goodness of fit indicators}
#' \item{\code{arg}}{some key parameters}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' ## The following code takes a long time to run:
#' econineq = sf::read_sf(system.file('extdata/econineq.gpkg',package = 'geocomplexity'))
#' gwr_geoc(formula = Gini ~ ., data = econineq,
#'          bw = "AIC", adaptive = TRUE)
#' }
gwr_geoc = \(formula, data, gcs = NULL, alpha = seq(0.05,1,0.05),
             bw = "RMSE", adaptive = TRUE, kernel = "gaussian"){
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)

  if (!inherits(data,'sf')){
    data = sf::st_as_sf(data)
  }

  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }

  yname = formula.vars[1]
  if (is.null(gcs)) {
    gcs = data %>%
      dplyr::select(-dplyr::any_of(yname)) %>%
      geocomplexity::geocd_vector(returnsf = FALSE) %>%
      as.matrix()
  }
  distm = sdsfun::sf_distance_matrix(data)
  geom = sf::st_geometry(data)
  data = sf::st_drop_geometry(data)
  xname = colnames(data)[-which(colnames(data) == yname)]
  y = data[,yname,drop = TRUE]
  xs = data %>%
    dplyr::select(dplyr::all_of(xname)) %>%
    as.matrix()

  bws = 0
  knns = 0
  criterion = "RMSE"
  if (adaptive) {
    if (is.character(bw)){
      criterion = bw
      knns = seq(floor(0.05*length(y)),ceiling(0.95*length(y)),by = 1)
    } else {
      bws = bw
    }
  } else {
    if (is.character(bw)){
      criterion = bw
      bws = SelectSortedBW(distm,0.1,0.5)
    } else {
      bws = bw
    }
  }

  optarg = SGWRSel(bws, knns, alpha, y, xs, distm,
                   adaptive, criterion, kernel)
  res = GeoCGWR(y,xs,gcs,distm,optarg[[1]],optarg[[2]],adaptive,optarg[[3]],kernel)
  res$SDF = purrr::map2_dfc(
                       res$SDF,
                       list(c("Intercept",xname),
                            paste0(c("Intercept",xname),"_SE"),
                            paste0(c("Intercept",xname),"_TV"),
                            c("Pred"),c("Residuals"),c("LocalR2")),
                       \(.mat,.name) {
                         colnames(.mat) = .name
                         return(tibble::as_tibble(.mat))
                       }) %>%
    sf::st_set_geometry(geom)
  res$arg = append(res$arg,list("criterion" = criterion))
  class(res) = 'gcgwrm'
  return(res)
}

#' @title print GCGWR result
#'
#' @param x Return by `gwr_geoc()`.
#' @param ... (optional) Other arguments.
#'
#' @return Formatted string output
#' @export
print.gcgwrm = \(x,...){
  coefdf = sf::st_drop_geometry(x$SDF) %>%
    dplyr::select(Intercept:Intercept_SE) %>%
    dplyr::select(-Intercept_SE)
  PrintGCGWRM(x,coefdf)
}

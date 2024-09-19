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
      geocd_vector(returnsf = FALSE) %>%
      as.matrix()
  }
  distm = sdsfun::sf_distance_matrix(data)
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
  return(res)
}

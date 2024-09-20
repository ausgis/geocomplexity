#include <RcppArmadillo.h>
#include "GWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat GeoCD_SWM(arma::mat gcs,
                    arma::mat wt,
                    std::string style){
  int n = wt.n_rows;
  arma::mat gc_wt = zeros(n,n);

  for (int i = 0; i < n; ++i) {
    arma::rowvec gcw =  zeros(n);
    for (int j = 0; j < n; ++j) {
      double gc = RowDiffAbsMean(gcs,i,j);
      gc = exp(-pow(gc,2));
      gc_wt(i,j) = gc;
    }
  }
  return gc_wt;
}

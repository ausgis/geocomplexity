#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCDependence(NumericVector x,
                                   NumericMatrix wt){
  NumericVector out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    double zi = x[i];
    NumericVector wtj = wt(i, _);
    IntegerVector j = rcpp_which(wtj != 0);
    int m = j.size();
    NumericVector zj = x[j];
    NumericVector wtj_z = multiply_vector(wtj[j],zj);
    double localf = -1.0 / m * zi * sum_nona(wtj_z);
    double surroundf = 0;
    for (int n = 0; n < m; ++n) {
      int kn = j[n];
      NumericVector wtk = wt(kn, _);
      IntegerVector k = rcpp_which(wtk != 0);
      k = intersect(j, k);
      NumericVector zk = x[k];
      NumericVector wtk_z = multiply_vector(wtj[k],zk);
      wtk_z = multiply_vector(wtk[k],wtk_z);
      double surroundv = mean_nona(wtk_z);
      surroundf += zj[n] * surroundv;
    }
    out[i] = localf + -1.0 / m * surroundf;
  }
  return out;
}

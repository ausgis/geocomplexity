#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCDependence(NumericVector x,
                                   IntegerMatrix wt){
  NumericVector out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    double zi = x[i];
    IntegerVector wtj = wt(i , _);
    IntegerVector j = rcpp_which(wtj != 0);
    int m = j.size();
    NumericVector zj = x[j];
    double localf = -1.0 / m * zi * sum_nona(zj);
    double surroundf = 0;
    for (int n = 0; n < m; ++n) {
      IntegerVector wtk = wt(j[n] , _);
      IntegerVector k = rcpp_which(wtk != 0);
      k = intersect(j, k);
      NumericVector zk = zj[k];
      double surroundv = mean_nona(zk);
      surroundf += zj[n] * surroundv;
    }
    out[i] = localf + -1.0 / m * surroundf;
  }
  return out;
}

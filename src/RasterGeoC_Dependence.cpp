#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double RasterGeoCDependenceOne(NumericVector x,
                               size_t ni = 9,
                               size_t nw = 9){
  NumericVector x1 = x[!is_na(x)];
  int m = x1.size() - 1;
  if (m <= 0) {
    m = 1;
  }
  int center_idx = (nw - 1) / 2;
  double zi = x[center_idx];
  NumericVector wt_i (nw,1);
  wt_i[center_idx] = 0;
  List wt_j = RasterGeoCNeighbors(sqrt(nw),sqrt(nw));
  NumericVector x_j = x[wt_i != 0];
  double localf = -1.0 / m * zi * sum_nona(x * wt_i);
  double surroundf = 0;
  for (int n = 0; n < wt_j.size(); ++n){
    IntegerVector j_index = wt_j[n];
    NumericVector xj = x[j_index];
    double surroundv = mean_nona(xj);
    double x_add = x_j[n] * surroundv;
    if (NumericVector::is_na(x_add)) {
      x_add = 0;
    }
    surroundf += x_add;
  }
  double res = localf + -1.0 / m * surroundf;
  return res;
}

// [[Rcpp::export]]
NumericVector RasterGeoCDependence(NumericVector x,
                                   size_t ni = 9,
                                   size_t nw = 9){
    NumericVector out(ni);
    for (size_t i = 0; i < ni; ++i) {
      NumericVector chunk(nw);
      for (size_t j = 0; j < nw; ++j) {
        chunk[j] = x[i * nw + j];
      }
      out[i] = RasterGeoCDependenceOne(chunk,ni,nw);
    }
    return out;
}

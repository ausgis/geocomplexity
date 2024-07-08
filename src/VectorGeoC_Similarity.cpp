#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCSimilarity(NumericVector x,
                                   NumericMatrix wt){
  NumericVector out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    double zi = x[i];
    NumericVector zs(x.size());
    for (int n = 0; n < x.size(); ++n) {
      zs[n] = cosine_similarity(zi,x[n]);
    }
    out[i] = spatial_variance(zs,wt);
  }
  return out;
}

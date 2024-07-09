#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCSimilarity(NumericMatrix xobs,
                                   NumericMatrix wt,
                                   int method){
  if (method != 1){
    NumericVector out(xobs.nrow());
    for (int i = 0; i < xobs.nrow(); ++i) {
      NumericVector zi = xobs(i,_);
      NumericVector zs(xobs.nrow());
      for (int n = 0; n < xobs.nrow(); ++n) {
        zs[n] = CosineSimilarity(zi,xobs(n,_));
      }
      out[i] = spatial_variance(zs,wt);
    }
    return out;
  } else {
    NumericVector out = GCS_Variance(xobs,wt);
    return out;
  }
}

#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCSimilarity(NumericMatrix xobs,
                                   NumericMatrix wt,
                                   int method){
  NumericVector out(xobs.nrow());
  if (method != 1){
    for (int i = 0; i < xobs.nrow(); ++i) {
      NumericVector zi = xobs(i,_);
      NumericVector zs(xobs.nrow());
      for (int n = 0; n < xobs.nrow(); ++n) {
        zs[n] = CosineSimilarity(zi,xobs(n,_));
      }
      out[i] = spatial_variance(zs,wt);
    }
  } else {

  }

  return out;
}

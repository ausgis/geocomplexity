#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCSimilarity(NumericMatrix xobs,
                                   NumericMatrix wt,
                                   int similarity,
                                   String method){
  if (similarity != 1){
    NumericVector out(xobs.nrow());
    for (int i = 0; i < xobs.nrow(); ++i) {
      NumericVector zi = xobs(i,_);
      NumericVector zs(xobs.nrow());
      for (int n = 0; n < xobs.nrow(); ++n) {
        zs[n] = CosineSimilarity(zi,xobs(n,_));
      }
      if (method == "spvar"){
        out[i] = spatial_variance(zs,wt);
      } else{
        out[i] = InforEntropy(zs);
      }
    }
    return out;
  } else {
    NumericVector out = GCS_Variance(xobs,wt,method);
    return out;
  }
}

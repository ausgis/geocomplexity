#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector VectorGeoCSSH(NumericMatrix xobs,
                            NumericMatrix wt,
                            String method){
  NumericVector out = SSH_Variance(xobs,wt,method);
  return out;
}

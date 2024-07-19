#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector RasterGeoCSSH(NumericMatrix x,
                            IntegerMatrix iw,
                            int w,String method){
  NumericMatrix wt (x.nrow(),x.nrow());
  IntegerVector w_sp = extract_window(iw,w);
  int ncell = pow(w,2);
  for (int n = 0; n < x.nrow(); ++n){
    IntegerVector wi = rcpp_seq(n*ncell,(n+1)*ncell-1);
    wi = w_sp[wi];
    wi = wi[!is_na(wi)];
    for (int ni = 0; ni < wi.size(); ++ni){
      wt(n,wi[ni]) = 1;
      wt(wi[ni],n) = 1;
    }
  }
  // Rcout << "wt:  "<< wt;
  NumericVector out = SSH_Variance(x,wt,method);
  return out;
}

#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector RasterGeoCSimilarity(NumericMatrix x,
                                   IntegerMatrix iw,
                                   int w){
  NumericVector out(x.nrow());
  IntegerVector w_sp = extract_window(iw,w);
  int ncell = pow(w,2);
  for (int n = 0; n < x.nrow(); ++n){
    IntegerVector wi = rcpp_seq(n*ncell,(n+1)*ncell-1);
    Rcout << "wi_1: " << wi << "\n";
    wi = w_sp[wi];
    wi = wi[!is_na(wi)];
    Rcout << "wi_2: " << wi << "\n";
    Rcout << "wi_2.size(): " << wi.size() << "\n";
    NumericVector zn = x(n,_);
    NumericVector resn(wi.size());
    for (int i = 0 ; i < wi.size(); ++i){
      NumericVector zi = x(wi[i],_);
      double zs = cosine_similarity(zn,zi);
      resn[i] = zs;
    }
    double res = 0;
    for(int i = 0; i < wi.size(); ++i) {
      for(int j = 0; j < wi.size(); ++j) {
        res += pow((resn[i]-resn[j]),2) / 2;
      }
    }
    out[n] = res / pow(wi.size(),2);
    Rcout << "out[n]: " << out[n] << "\n";
  }
  return out;
}

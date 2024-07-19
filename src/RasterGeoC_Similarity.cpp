#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector RasterGeoCSimilarity(NumericMatrix x, IntegerMatrix iw,
                                   int w,int similarity,String method){
  if (similarity != 1){
    NumericVector out(x.nrow());
    IntegerVector w_sp = extract_window(iw,w);
    int ncell = pow(w,2);
    for (int n = 0; n < x.nrow(); ++n){
      IntegerVector wi = rcpp_seq(n*ncell,(n+1)*ncell-1);
      // Rcout << "wi_1: " << wi << "\n";
      wi = w_sp[wi];
      // Rcout << "wi_2: " << wi << "\n";
      wi = wi[!is_na(wi)];
      // Rcout << "wi_3: " << wi << "\n";
      // Rcout << "wi_3.size(): " << wi.size() << "\n";
      NumericVector zn = x(n,_);
      NumericVector resn(wi.size());
      for (int i = 0 ; i < wi.size(); ++i){
        NumericVector zi = x(wi[i],_);
        double zs = CosineSimilarity(zn,zi);
        resn[i] = zs;
      }
      if (method == "spvar"){
        double res = 0;
        for(int i = 0; i < wi.size(); ++i) {
          for(int j = 0; j < wi.size(); ++j) {
            res += pow((resn[i]-resn[j]),2) / 2;
          }
        }
        out[n] = res / pow(wi.size(),2);
      } else {
        out[n] = InforEntropy(resn);
      }
    }
    return out;
  } else {
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
    NumericVector out = GCS_Variance(x,wt,method);
    return out;
  }
}

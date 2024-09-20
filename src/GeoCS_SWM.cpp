#include <Rcpp.h>
#include "GeoC_Helper.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix GeoCS_SWM(NumericVector x,
                        NumericMatrix wt,
                        String style){
  NumericMatrix result(x.size(), x.size());
  x = NormalizeVector(x,1,std::exp(1.0));
  for (int i = 0; i < x.size(); ++i) {
    for (int j = 0; j < x.size(); ++j) {
      double wij = wt(i,j);
      if (wij != 0) {
        double z = x[i] / x[j];
        result(i,j) = std::log(z);
      } else {
        result(i,j) = 0;
      }
    }
  }

  if (style == "B"){
    result = NormalizeMatRow(result,0,1);
  }
  if (style == "W"){
    result = MatRowStandardize(result);
  }
  if (style == "C"){
    result = MatGlobalStandardize(result);
  }

  return result;
}

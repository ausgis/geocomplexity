#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

NumericMatrix GeoCSWT(NumericVector x,
                      NumericMatrix wt){
  NumericMatrix result(x.size(), x.size());
  for (int i = 0; i < x.size(); ++i) {
    for (int j = 0; j < x.size(); ++j) {
      double wij = wt(i,j);
      if (wij != 0) {
        result(i,j) = x[i] / x[j];
      } else {
        result(i,j) = 0;
      }
    }
  }
  return result;
}

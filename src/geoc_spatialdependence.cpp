#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double geoc_spatialdependence(NumericVector x){
    NumericVector wt_i = NumericVector::create(1,1,1,1,0,1,1,1,1);
    List wt_j = List::create(IntegerVector::create(1,3),
                             IntegerVector::create(0,2,3,5),
                             IntegerVector::create(1,5),
                             IntegerVector::create(0,1,6,7),
                             IntegerVector::create(1,2,7,8),
                             IntegerVector::create(3,7),
                             IntegerVector::create(3,5,6,8),
                             IntegerVector::create(5,7));
    NumericVector x_j = x[wt_i != 0];
    double localf = -1.0 / 8 * x[4] * Rcpp::sum(x * wt_i);
    double surroundf = 0;
    for(int n = 0; n < wt_j.size(); ++n){
      IntegerVector j_index = wt_j[n];
      NumericVector xj = x[j_index];
      double surroundv = Rcpp::mean(xj);
      surroundf += x_j[n] * surroundv;
    }
    double res = localf + -1.0 / 8 * surroundf;
    return res;
}

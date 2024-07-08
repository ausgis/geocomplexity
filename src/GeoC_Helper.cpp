#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double sum_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  return sum(x1);
}

double mean_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  return mean(x1);
}

IntegerVector rcpp_which(LogicalVector x){
  IntegerVector v = seq(0,x.size()-1);
  return v[x];
}

NumericVector multiply_vector(IntegerVector intVec, NumericVector numVec) {
  int n = intVec.size();
  NumericVector result(n);

  for(int i = 0; i < n; ++i) {
    result[i] = intVec[i] * numVec[i];
  }

  return result;
}

List remove_index(List lst, int idx) {
  int n = lst.size();
  if (idx < 0 || idx >= n) {
    stop("Index out of bounds");
  }
  List result(n - 1);
  int result_idx = 0;
  for (int i = 0; i < n; ++i) {
    if (i != idx) {
      result[result_idx++] = lst[i];
    }
  }
  return result;
}

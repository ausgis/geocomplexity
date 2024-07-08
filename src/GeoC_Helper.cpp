#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double sum_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  double xres = sum(x1);
  // if (NumericVector::is_na(xres)) {
  //   xres = 0;
  // }
  return xres;
}

double mean_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  double xres = mean(x1);
  // if (NumericVector::is_na(xres)) {
  //   xres = 0;
  // }
  return xres;
}

double sd_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  double xres = sd(x1);
  // if (NumericVector::is_na(xres)) {
  //   xres = 0;
  // }
  return xres;
}

double cosine_similarity(NumericVector A, NumericVector B) {
  int n = A.size();
  double dot_product = 0.0;
  double norm_A = 0.0;
  double norm_B = 0.0;

  for(int i = 0; i < n; i++) {
    dot_product += A[i] * B[i];
    norm_A += A[i] * A[i];
    norm_B += B[i] * B[i];
  }

  norm_A = sqrt(norm_A);
  norm_B = sqrt(norm_B);

  return dot_product / (norm_A * norm_B);
}

IntegerVector rcpp_which(LogicalVector x){
  IntegerVector v = seq(0,x.size()-1);
  return v[x];
}

NumericVector multiply_vector(NumericVector numvec1, NumericVector numvec2) {
  int n = numvec1.size();
  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = numvec1[i] * numvec2[i];
  }
  return result;
}

double matrix_sum(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  double out = 0;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      out += mat(i, j);
    }
  }
  return out;
}

// [[Rcpp::export]]
double spatial_variance(NumericVector x, NumericMatrix wt) {
  int n = x.size();
  double out = 0;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      double w = wt(i, j);
      out += w * pow((x[i]-x[j]),2) / 2;
    }
  }
  out = out / matrix_sum(wt);
  return out;
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

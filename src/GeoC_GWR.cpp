#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double gaussian_kernel(double dist, double bw) {
  return (1.0 / (sqrt(2.0 * M_PI) * bw)) * exp(-0.5 * pow(dist / bw, 2));
}

// [[Rcpp::export]]
mat GeoCGWRFit(arma::vec y, arma::mat X, arma::vec gcs, arma::mat Cdist,
               double bw, String kernel = "gaussian", double alpha = 0.5) {
  int n = X.n_rows;
  int k = X.n_cols;
  mat X_with_intercept = join_horiz(ones<mat>(n, 1), X);
  mat betas = zeros(n, k + 1);

  for (int i = 0; i < n; ++i) {
    mat gc_wt = zeros(n);
    vec dist_wt = zeros(n);

    // Calculate weight matrix
    for (int j = 0; j < n; ++j) {
      double dist = Cdist(i,j);
      double gc = gcs(j) / gcs(i);
      gc_wt(j) = gc;
      dist_wt(j) = gaussian_kernel(dist, bw);
    }
    vec wt = alpha * gc_wt + (1 - alpha) * dist_wt;
    mat W = diagmat(wt);
    // Weighted Least Squares
    mat XtWX = X_with_intercept.t() * W * X_with_intercept;
    vec XtWy = X_with_intercept.t() * W * y;

    // Solve local regression coefficient
    betas.row(i) = solve(XtWX, XtWy).t();
  }
  return betas;
}


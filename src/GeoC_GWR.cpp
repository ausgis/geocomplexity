#include <RcppArmadillo.h>
#include "GeoCGWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat GeoCGWRFit(arma::vec y, arma::mat X, arma::vec gcs, arma::mat Cdist,
                     double bw, std::string kernel = "boxcar", double alpha = 0.5) {
  int n = X.n_rows;
  int k = X.n_cols;
  arma::mat X_with_intercept = arma::join_horiz(ones<mat>(n, 1), X);
  arma::mat betas = arma::zeros(n, k + 1);

  for (int i = 0; i < n; ++i) {
    arma::mat gc_wt = arma::zeros(n);
    arma::vec dist_wt = arma::zeros(n);

    // Calculate weight matrix
    for (int j = 0; j < n; ++j) {
      double dist = Cdist(i,j);
      double gc = gcs(j) / gcs(i);
      gc_wt(j) = gc;

      if (kernel == "gaussian") {
        dist_wt(j) = gaussian_kernel(dist, bw);
      } else if (kernel == "exponential") {
        dist_wt(j) = exponential_kernel(dist, bw);
      } else if (kernel == "bisquare") {
        dist_wt(j) = bisquare_kernel(dist, bw);
      } else if (kernel == "triangular") {
        dist_wt(j) = triangular_kernel(dist, bw);
      }  else if (kernel == "boxcar") {
        dist_wt(j) = boxcar_kernel(dist, bw);
      }

    }
    arma::vec wt = alpha * gc_wt + (1 - alpha) * dist_wt;
    arma::mat W = arma::diagmat(wt);
    // Weighted Least Squares
    arma::mat XtWX = X_with_intercept.t() * W * X_with_intercept;
    arma::vec XtWy = X_with_intercept.t() * W * y;

    // Solve local regression coefficient
    betas.row(i) = arma::solve(XtWX, XtWy).t();
  }

  return betas;
}


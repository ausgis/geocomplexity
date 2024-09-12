#include <RcppArmadillo.h>
#include "GeoCGWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List GeoCGWRFit(arma::vec y, arma::mat X,
                      arma::vec gcs, arma::mat Cdist,
                      double bw, std::string kernel = "gaussian") {
  int n = X.n_rows;
  int k = X.n_cols;
  arma::mat X_with_intercept = arma::join_horiz(ones<mat>(n, 1), X);
  arma::mat betas = arma::zeros(n, k + 1);
  arma::mat se_betas = zeros(n, k + 1);
  arma::mat hat_matrix = zeros(n,n);
  arma::vec residuals = zeros(n);
  arma::vec yhat = zeros(n);

  for (int i = 0; i < n; ++i) {
    arma::vec gc_wt = arma::zeros(n);
    arma::vec dist_wt = arma::zeros(n);

    // Calculate Weight Matrix
    for (int j = 0; j < n; ++j) {
      double gc = gcs(j) / gcs(i);
      gc_wt(j) = gc;

      double dist = Cdist(i,j);
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
    arma::vec wt = gc_wt % dist_wt;
    arma::mat W = arma::diagmat(wt);
    // Weighted Least Squares
    arma::mat XtWX = X_with_intercept.t() * W * X_with_intercept;
    arma::vec XtWy = X_with_intercept.t() * W * y;

    // Regularization to avoid singular matrix
    arma::mat XtWX_reg = XtWX + 1e-5 * arma::eye(XtWX.n_rows, XtWX.n_cols);

    // Solve Local Regression Coefficient
    arma::vec beta_i = arma::solve(XtWX_reg, XtWy);
    betas.row(i) = beta_i.t();

    // Calculate Residuals: y_i - X_i * beta_i
    double y_hat_i = arma::as_scalar(X_with_intercept.row(i) * beta_i);
    residuals(i) = y(i) - y_hat_i;
    yhat(i) = y_hat_i;

    // Calculate the Diagonal Elements of the Cap Hat Matrix
    arma::mat XtWX_inv = arma::inv(XtWX);
    arma::rowvec hat_row = X_with_intercept.row(i) * XtWX_inv * X_with_intercept.t() * W;
    hat_matrix.row(i) = hat_row;

    // Standard Error of Calculated Coefficient
    double sigma2_i = arma::as_scalar(sum(pow(residuals(i), 2)) / (n - k - 1));
    arma::vec se_beta_i = sqrt(sigma2_i * arma::diagvec(XtWX_inv));
    se_betas.row(i) = se_beta_i.t();
  }

  arma::mat t_values = betas / se_betas;
  double rss = arma::as_scalar(arma::sum(arma::pow(residuals, 2)));
  double tss = arma::as_scalar(arma::sum(arma::pow(y - arma::mean(y), 2)));
  double r2 = 1 - (rss / tss);
  double adjr2 = 1 - (1 - r2) * (n - 1) / (n - k - 1);
  double rmse = sqrt(rss / n);
  double aic = n * log(rss / n) + 2 * k;
  double sigma_hat = rss / (n - (2*arma::trace(hat_matrix)-arma::trace(hat_matrix.t()*hat_matrix)));
  double aicc = 2*n*log(sigma_hat) + n*log(2*M_PI) + n*(n+arma::trace(hat_matrix)/(n-2-arma::trace(hat_matrix)));

  return Rcpp::List::create(
    Named("Coefficient") = betas,
    Named("t_values") = t_values,
    Named("SE_Coefficient") = se_betas,
    Named("yhat") = yhat,
    Named("Residuals") = residuals,
    Named("RSS") = rss,
    Named("R2") = r2,
    Named("R2_Adj") = adjr2,
    Named("RMSE") = rmse,
    Named("AIC") = aic,
    Named("AICc") = aicc
  );
}



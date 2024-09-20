#include <RcppArmadillo.h>
#include "GWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List GeoCGWRFit(arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist,
                      double bw = 0, double knn = 0, bool adaptive = false,
                      double alpha = 0.05, std::string kernel = "gaussian") {
  int n = X.n_rows;
  int k = X.n_cols;
  arma::vec bw_vec;
  arma::mat X_sd = StandardizeMatColumns(X);
  arma::mat X_with_intercept = arma::join_horiz(ones<mat>(n, 1), X);
  arma::mat betas = arma::zeros(n, k + 1);
  arma::mat se_betas = zeros(n, k + 1);
  arma::mat hat_matrix = zeros(n,n);
  arma::vec residuals = zeros(n);
  arma::vec yhat = zeros(n);
  arma::vec local_r2_vec = arma::zeros(n);

  if (adaptive) {
    bw_vec = GenAdaptiveKNNBW(Gdist,knn);
  } else {
    bw_vec = Double4Vec(bw);
  }

  for (int i = 0; i < n; ++i) {
    arma::vec gc_wt = arma::zeros(n);
    arma::vec dist_wt = arma::zeros(n);

    // Determine bandwidth for the current point
    double current_bw = adaptive ? bw_vec(i) : bw_vec(0);

    // Calculate Weight Matrix
    for (int j = 0; j < n; ++j) {
      double gc = RowDiffAbsMeanWeight(X_sd,arma::trans(gcs.row(j)),i,j);
      gc_wt(j) = exp(-pow(gc,2));

      double dist = Gdist(i,j);
      if (kernel == "gaussian") {
        dist_wt(j) = gaussian_kernel(dist, current_bw);
      } else if (kernel == "exponential") {
        dist_wt(j) = exponential_kernel(dist, current_bw);
      } else if (kernel == "bisquare") {
        dist_wt(j) = bisquare_kernel(dist, current_bw);
      } else if (kernel == "triangular") {
        dist_wt(j) = triangular_kernel(dist, current_bw);
      }  else if (kernel == "boxcar") {
        dist_wt(j) = boxcar_kernel(dist, current_bw);
      }

    }
    arma::vec wt = alpha * gc_wt + (1 - alpha) * dist_wt;
    // wt = Normalize4Interval(wt,0,1);
    arma::mat W = arma::diagmat(wt);

    // Weighted Least Squares
    arma::mat XtWX = X_with_intercept.t() * W * X_with_intercept;
    arma::vec XtWy = X_with_intercept.t() * W * y;

    // Regularization to avoid singular matrix
    arma::mat XtWX_reg = XtWX + 1e-6 * arma::eye(XtWX.n_rows, XtWX.n_cols);

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

    // Local R-square calculation
    local_r2_vec(i) = LocalR2(y, X_with_intercept * beta_i, wt);

    // Standard Error of Calculated Coefficient
    double sigma2_i = arma::as_scalar(sum(pow(residuals(i), 2)) / (n - k - 1));
    arma::vec se_beta_i = sqrt(sigma2_i * arma::diagvec(XtWX_inv));
    se_betas.row(i) = se_beta_i.t();
  }

  // Compute additional metrics
  // refer to https://github.com/rsbivand/spgwr/blob/main/R/gwr.R
  arma::mat t_values = betas / se_betas;
  double rss = arma::as_scalar(arma::sum(arma::pow(residuals, 2)));
  double tss = arma::as_scalar(arma::sum(arma::pow(y - arma::mean(y), 2)));
  double r2 = 1 - (rss / tss);
  double adjr2 = 1 - (1 - r2) * (n - 1) / (n - k - 1);
  double rmse = sqrt(rss / n);
  double v1 = arma::trace(hat_matrix);
  arma::mat B2 = hat_matrix.t() * hat_matrix;
  double v2 = arma::trace(B2);
  // effective n.p. is 2*v1 - v2
  double enp = 2*v1 - v2;
  // effective d.f. is n - 2*v1 + v2
  double edf = n - 2*v1 + v2;
  arma::mat B1 = (DiagMatrix(n) - hat_matrix).t() * (DiagMatrix(n) - hat_matrix);
  double delta1 = arma::trace(B1);
  // double sigma2 = rss / delta1;
  // double odelta2 = sum(arma::pow(arma::diagvec(B1), 2));
  double delta2 = sum(arma::diagvec(B1*B1));
  double nu1 = arma::trace(B2);
  double sigma2_b = rss / n;
  // double aicbb = 2*n*log(sqrt(sigma2_b)) + n*log(2*datum::pi) + (n * ((n + v1) / (n - 2 - v1)));
  double aichb = 2*n*log(sqrt(sigma2_b)) + n*log(2*datum::pi) + n + v1;
  double aiccb = 2*n*log(sqrt(sigma2_b)) + n*log(2*datum::pi) + n * ((delta1/delta2)*(n + nu1))/((pow(delta1,2)/delta2)-2);

  return Rcpp::List::create(
    Named("Coefficient") = betas,
    Named("SE_Coefficient") = se_betas,
    Named("t_values") = t_values,
    Named("Pred") = yhat,
    Named("Residuals") = residuals,
    Named("RSS") = rss,
    Named("ENP") = enp,
    Named("EDF") = edf,
    Named("R2") = r2,
    Named("R2_Adj") = adjr2,
    Named("RMSE") = rmse,
    // Named("AICb") = aicbb,
    Named("AIC") = aichb,
    Named("AICc") = aiccb,
    Rcpp::Named("LocalR2") = local_r2_vec
  );
}

// [[Rcpp::export]]
Rcpp::List GeoCGWR(arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist,
                   double bw = 0, double knn = 0, bool adaptive = false,
                   double alpha = 0.01, std::string kernel = "gaussian") {
  Rcpp::List g = GeoCGWRFit(y,X,gcs,Gdist,bw,knn,adaptive,alpha,kernel);
  return Rcpp::List::create(
    Rcpp::Named("SDF") = Rcpp::List::create(
      Rcpp::Named("Coefficient") = g["Coefficient"],
      Rcpp::Named("SE_Coefficient") = g["SE_Coefficient"],
      Rcpp::Named("t_values") = g["t_values"],
      Rcpp::Named("Pred") = g["Pred"],
      Rcpp::Named("Residuals") = g["Residuals"],
      Rcpp::Named("LocalR2") = g["LocalR2"]),
    Rcpp::Named("diagnostic") = Rcpp::List::create(
          Rcpp::Named("RSS") = g["RSS"],
          Rcpp::Named("ENP") = g["ENP"],
          Rcpp::Named("EDF") = g["EDF"],
          Rcpp::Named("R2") = g["R2"],
          Rcpp::Named("R2_Adj") = g["R2_Adj"],
          Rcpp::Named("RMSE") = g["RMSE"],
          Rcpp::Named("AIC") = g["AIC"],
          Rcpp::Named("AICc") = g["AICc"]),
     Rcpp::Named("arg") = Rcpp::List::create(
       Rcpp::Named("adaptive") = adaptive,
       Rcpp::Named("bw") = bw,
       Rcpp::Named("knn") = knn,
       Rcpp::Named("alpha") = alpha)
  );
}

#include <RcppArmadillo.h>
#include "GeoCGWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List GeoCSGWRFit(arma::vec y, arma::mat X, arma::vec gcs, arma::mat Gdist,
                       double bw  = 0, double knn = 0, bool adaptive = false,
                       double alpha = 0.5, std::string kernel = "gaussian") {
  int n = X.n_rows;
  int k = X.n_cols;
  arma::vec bw_vec;
  arma::mat X_with_intercept = arma::join_horiz(ones<mat>(n, 1), X);
  arma::mat betas = arma::zeros(n, k + 1);
  arma::mat se_betas = zeros(n, k + 1);
  arma::mat hat_matrix = zeros(n,n);
  arma::vec residuals = zeros(n);
  arma::vec yhat = zeros(n);

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
      double gc = gcs[i] - gcs[j];
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
    wt = Normalize4Interval(wt,0,1);
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
    Named("AICc") = aiccb
  );
}

// [[Rcpp::export]]
Rcpp::List GeoCGWRSel(arma::vec bandwidth, arma::vec knns, arma::vec alpha,
                      arma::vec y, arma::mat X, arma::vec gcs, arma::mat Cdist,
                      bool adaptive = false, std::string kernel = "gaussian") {
  if (adaptive) {
    int n = knns.n_elem;
    int k = alpha.n_elem;
    double AIC = std::numeric_limits<double>::max();
    double opt_knn = 0;
    double opt_alpha = 0;

    for (int i = 0; i < n; ++i) {
      double knn = knns(i);
      for (int j = 0; j < k; ++j) {
        double alpha_sel = alpha(j);
        List GeoCGWRResult = GeoCGWRFit(y,X,gcs,Cdist,0,knn,true,alpha_sel,kernel);
        double AICSel = GeoCGWRResult["AICc"];
        if (AICSel < AIC) {
          AIC = AICSel;
          opt_knn = knn;
          opt_alpha = alpha_sel;
        }
      }
    }
    return Rcpp::List::create(
      Named("bw") = 0,
      Named("knn") = opt_knn,
      Named("alpha") = opt_alpha
    );
  } else {
    int n = bandwidth.n_elem;
    int k = alpha.n_elem;
    double AIC = std::numeric_limits<double>::max();
    double opt_bw = 0;
    double opt_alpha = 0;

    for (int i = 0; i < n; ++i) {
      double bw = bandwidth(i);
      for (int j = 0; j < k; ++j) {
        double alpha_sel = alpha(j);
        List GeoCGWRResult = GeoCGWRFit(y,X,gcs,Cdist,bw,0,false,alpha_sel,kernel);
        double AICSel = GeoCGWRResult["AICc"];
        if (AICSel < AIC) {
          AIC = AICSel;
          opt_bw = bw;
          opt_alpha = alpha_sel;
        }
      }
    }
    return Rcpp::List::create(
      Named("bw") = opt_bw,
      Named("knn") = 0,
      Named("alpha") = opt_alpha
    );
  }
}

// [[Rcpp::export]]
Rcpp::List GeoCGWR(arma::vec y, arma::mat X, arma::vec gcs, arma::mat Cdist,
                   SEXP bw, arma::vec alpha, bool adaptive = false,
                   std::string kernel = "gaussian") {
  arma::vec knns;
  arma::vec bws;
  double MaxD = MaxInMatrix(Cdist);
  double MinD = MinInMatrix(Cdist);
  if (TYPEOF(bw) == STRSXP) {
    if (adaptive) {
      knns = ArmaSeq(3,15,1);
      bws = Double4Vec(0);
    } else {
      knns = Double4Vec(0);
      bws = ArmaSeq(MinD,MaxD,13);
    }
  } else if (TYPEOF(bw) == REALSXP) {
    NumericVector numericInput(bw);
    arma::vec v(numericInput.size());
    for (int i = 0; i < numericInput.size(); ++i) {
      v[i] = numericInput[i];
    }
    bws = v;
  } else {
    stop("Unsupported input type.");
  }

  Rcpp::List res = GeoCGWRSel(bws,knns,alpha,y,X,gcs,Cdist,adaptive,kernel);
  double optbw = res[0];
  double optknn = res[1];
  double optalpha = res[2];
  res = GeoCGWRFit(y,X,gcs,Cdist,optbw,optknn,adaptive,optalpha,kernel);
  return res;
}
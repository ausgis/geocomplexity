#include <RcppArmadillo.h>
#include "GeoCGWR_Helper.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List GeoCGWRFit(arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist,
                      double bwc = 0, double bwg = 0, double knn = 0,
                      bool adaptive = false, double alpha = 0.5,
                      std::string kernel = "gaussian") {
  int n = X.n_rows;
  int k = X.n_cols;
  arma::vec bwc_vec;
  arma::vec bwg_vec;
  arma::mat X_with_intercept = arma::join_horiz(ones<mat>(n, 1), X);
  arma::mat betas = arma::zeros(n, k + 1);
  arma::mat se_betas = zeros(n, k + 1);
  arma::mat hat_matrix = zeros(n,n);
  arma::vec residuals = zeros(n);
  arma::vec yhat = zeros(n);
  for (int i = 0; i < k; ++i) {
    gcs.col(i) = Normalize4Interval(gcs.col(i),1,100000);
  }
  arma::mat Cdist = EucdistM(gcs);

  if (adaptive) {
    bwc_vec = GenAdaptiveKNNBW(Cdist,knn);
    bwg_vec = GenAdaptiveKNNBW(Gdist,knn);
  } else {
    bwc_vec = Double4Vec(bwc);
    bwg_vec = Double4Vec(bwg);
  }

    for (int i = 0; i < n; ++i) {
      arma::vec gc_wt = arma::zeros(n);
      arma::vec dist_wt = arma::zeros(n);

      // Determine bandwidth for the current point
      double current_bw1 = adaptive ? bwc_vec(i) : bwc_vec(0);
      double current_bw2 = adaptive ? bwg_vec(i) : bwg_vec(0);

      // Calculate Weight Matrix
      for (int j = 0; j < n; ++j) {
        double dist1 = Cdist(i,j);
        double dist2 = Gdist(i,j);
        if (kernel == "gaussian") {
          gc_wt(j) = gaussian_kernel(dist1, current_bw1);
          dist_wt(j) = gaussian_kernel(dist2, current_bw2);
        } else if (kernel == "exponential") {
          gc_wt(j) = exponential_kernel(dist1, current_bw1);
          dist_wt(j) = exponential_kernel(dist2, current_bw2);
        } else if (kernel == "bisquare") {
          gc_wt(j) = bisquare_kernel(dist1, current_bw1);
          dist_wt(j) = bisquare_kernel(dist2, current_bw2);
        } else if (kernel == "triangular") {
          gc_wt(j) = triangular_kernel(dist1, current_bw1);
          dist_wt(j) = triangular_kernel(dist2, current_bw2);
        } else if (kernel == "boxcar") {
          gc_wt(j) = boxcar_kernel(dist2, current_bw1);
          dist_wt(j) = boxcar_kernel(dist2, current_bw2);
        }

      }
      arma::vec wt = alpha * gc_wt + (1 - alpha) * dist_wt;
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
Rcpp::List GeoCGWRSel(arma::vec bwcs, arma::vec bwgs, arma::vec knns, arma::vec alpha,
                      arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist,
                      bool adaptive = false, std::string criterion = "RMSE",
                      std::string kernel = "gaussian") {
    if (adaptive) {
      int n = knns.n_elem;
      int k = alpha.n_elem;
      double measures = std::numeric_limits<double>::max();
      double opt_knn = 0;
      double opt_alpha = 0;

      for (int i = 0; i < n; ++i) {
        double knn = knns(i);
        for (int j = 0; j < k; ++j) {
          double alpha_sel = alpha(j);
          List GeoCGWRResult = GeoCGWRFit(y,X,gcs,Gdist,0,0,knn,true,alpha_sel,kernel);
          double MeaSel = GeoCGWRResult[criterion];
          if (MeaSel < measures) {
            measures = MeaSel;
            opt_knn = knn;
            opt_alpha = alpha_sel;
          }
        }
      }
      return Rcpp::List::create(
        Named("bwc") = 0,
        Named("bwg") = 0,
        Named("knn") = opt_knn,
        Named("alpha") = opt_alpha
      );
    } else {
      int n = bwcs.n_elem;
      int k = alpha.n_elem;
      double measures = std::numeric_limits<double>::max();
      double opt_bwc = 0;
      double opt_bwg = 0;
      double opt_alpha = 0;

      for (int i = 0; i < n; ++i) {
        double bwc = bwcs(i);
        double bwg = bwgs(i);
        for (int j = 0; j < k; ++j) {
          double alpha_sel = alpha(j);
          List GeoCGWRResult = GeoCGWRFit(y,X,gcs,Gdist,bwc,bwg,0,false,alpha_sel,kernel);
          double MeaSel = GeoCGWRResult[criterion];
          if (MeaSel < measures) {
            measures = MeaSel;
            opt_bwc = bwc;
            opt_bwg = bwg;
            opt_alpha = alpha_sel;
          }
        }
      }
      return Rcpp::List::create(
        Named("bwc") = opt_bwc,
        Named("bwg") = opt_bwg,
        Named("knn") = 0,
        Named("alpha") = opt_alpha
      );
    }
  }

// [[Rcpp::export]]
Rcpp::List GeoCGWR(arma::vec y, arma::mat X, arma::mat gcs,
                   arma::mat Gdist, SEXP bwc, SEXP bwg,
                   arma::vec alpha, bool adaptive = false,
                   std::string kernel = "gaussian") {
  arma::vec knns;
  arma::vec bwcs;
  arma::vec bwgs;
  arma::mat Cdist = EucdistM(gcs);
  std::string criterion = "RMSE";
  int sample_n = static_cast<int>(y.n_elem / 10);
  double MaxCD = MaxInMatrix(Cdist);
  double MinCD = MinInMatrix(Cdist);
  double MaxGD = MaxInMatrix(Gdist);
  double MinGD = MinInMatrix(Gdist);
  if (TYPEOF(bwc) == STRSXP) {
    if (adaptive) {
      knns = ArmaSeq(3,sample_n,1);
      bwcs = Double4Vec(0);
    } else {
      knns = Double4Vec(0);
      bwcs = arma::linspace(MinCD,MaxCD/3,sample_n);
    }
  } else if (TYPEOF(bwc) == REALSXP) {
    double v = Rcpp::as<double>(bwc);
    bwcs = Double4Vec(v);
  } else {
    stop("Unsupported input type.");
  }

  if (TYPEOF(bwg) == STRSXP) {
    std::string criterion = Rcpp::as<std::string>(bwc);
    if (adaptive) {
      knns = ArmaSeq(3,sample_n,1);
      bwgs = Double4Vec(0);
    } else {
      knns = Double4Vec(0);
      bwgs = arma::linspace(MinGD,MaxGD/3,sample_n);
    }
  } else if (TYPEOF(bwg) == REALSXP) {
    double v = Rcpp::as<double>(bwg);
    bwgs = Double4Vec(v);
  } else {
    stop("Unsupported input type.");
  }

  Rcpp::List res = GeoCGWRSel(bwcs,bwgs,knns,alpha,y,X,gcs,Gdist,adaptive,criterion,kernel);
  double optbwc = res[0];
  double optbwg = res[1];
  double optknn = res[2];
  double optalpha = res[3];
  res = GeoCGWRFit(y,X,gcs,Gdist,optbwc,optbwg,optknn,adaptive,optalpha,kernel);
  return res;
}

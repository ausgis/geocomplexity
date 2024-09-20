#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Gaussian Kernel
double gaussian_kernel(double dist, double bw) {
  return exp(-pow(dist / bw, 2));
}

// double gaussian_kernel(double dist, double bw) {
//   return exp(-0.5 * pow(dist / bw, 2));
// }

// Exponential Kernel
double exponential_kernel(double dist, double bw) {
  return exp(-dist / bw);
}

// Bisquare Kernel
double bisquare_kernel(double dist, double bw) {
  if (dist < bw) {
    return pow(1 - pow(dist / bw, 2), 2);
  } else {
    return 0.0;
  }
}

// Triangular Kernel
double triangular_kernel(double dist, double bw) {
  if (dist < bw) {
    return pow(1 - pow(dist / bw, 3), 3);
  } else {
    return 0.0;
  }
}

// Boxcar Kernel
double boxcar_kernel(double dist, double bw) {
  if (dist < bw) {
    return 1.0;
  } else {
    return 0.0;
  }
}

arma::vec Normalize4Interval(const arma::vec& v, double a, double b) {
  double min_v = arma::min(v);
  double max_v = arma::max(v);

  // Prevent division by zero
  if (max_v == min_v) {
    return arma::vec(v.n_elem, arma::fill::ones) * (a + b) / 2;
  }

  // Normalize to [0, 1]
  arma::vec v_normalized = (v - min_v) / (max_v - min_v);

  // Normalize to [a, b]
  arma::vec v_scaled = a + v_normalized * (b - a);

  return v_scaled;
}

arma::mat DiagMatrix(int n) {
  // Create a vector of ones with length n
  arma::vec diag_elements = arma::ones<arma::vec>(n);

  // Create a diagonal matrix using the vector
  arma::mat diag_mat = arma::diagmat(diag_elements);

  return diag_mat;
}

arma::vec ArmaSeq(double from, double to, double by = 1, int length_out = -1) {
  if (length_out > 0) {
    return arma::linspace(from, to, length_out);
  } else {
    if (by == 0) {
      Rcpp::stop("The `by` parameter cannot be 0.");
    }

    int n = static_cast<int>((to - from) / by) + 1;
    arma::vec result(n);
    for (int i = 0; i < n; i++) {
      result[i] = from + i * by;
    }
    return result;
  }
}

arma::vec Double4Vec(double x) {
  arma::vec v(1);
  v[0] = x;
  return v;
}

arma::vec GenAdaptiveKNNBW(const arma::mat& D, double k) {
  int n = D.n_rows;
  arma::vec bandwidths(n);

  for (int i = 0; i < n; ++i) {
    arma::rowvec distances = D.row(i);
    arma::uvec sorted_indices = sort_index(distances);
    double bandwidth = distances(sorted_indices(k));
    bandwidths(i) = bandwidth;
  }

  return bandwidths;
}

double MaxInMatrix(const arma::mat& mat) {
  double maxVal = arma::max(arma::max(mat));
  return maxVal;
}

double MinInMatrix(const arma::mat& mat) {
  double maxVal = arma::min(arma::min(mat));
  return maxVal;
}

arma::mat EucdistM(const arma::mat& X) {
  int n = X.n_rows;
  arma::mat dist_matrix(n, n, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = arma::norm(X.row(i) - X.row(j), 2);
      dist_matrix(i, j) = dist;
      dist_matrix(j, i) = dist;
    }
  }
  return dist_matrix;
}

double RowDiffAbsMean(const arma::mat& X, int i, int j) {
  arma::rowvec row_i = X.row(i);
  arma::rowvec row_j = X.row(j);
  arma::rowvec diff = arma::abs(row_i - row_j);
  return arma::mean(diff);
}

double MeanWeight(const arma::vec& x, const arma::vec& w) {
  double weighted_sum = arma::dot(x, w);
  double sum_weights = arma::sum(w);
  return weighted_sum / sum_weights;
}

double RowDiffAbsMeanWeight(const arma::mat& X,
                            const arma::vec& W,
                            int i, int j) {
  arma::rowvec row_i = X.row(i);
  arma::rowvec row_j = X.row(j);
  arma::rowvec abs_diff = arma::abs(row_i - row_j);
  double weighted_mean = arma::dot(abs_diff, W) / arma::sum(W);
  return weighted_mean;
}

arma::mat StandardizeMatColumns(const arma::mat& X) {
  arma::mat X_std = X;
  for (arma::uword i = 0; i < X.n_cols; ++i) {
    double mean_col = arma::mean(X.col(i));
    double std_col = arma::stddev(X.col(i));

    if (std_col != 0) {
      X_std.col(i) = (X.col(i) - mean_col) / std_col;
    } else {
      X_std.col(i).fill(0);
    }
  }

  return X_std;
}

Rcpp::NumericVector ArmaVec4RcppNumericVector(const arma::vec& arma_vec) {
  Rcpp::NumericVector numeric_vec(arma_vec.n_elem);
  std::copy(arma_vec.begin(), arma_vec.end(), numeric_vec.begin());

  return numeric_vec;
}

double LocalR2(const arma::vec& y, const arma::vec& yhat, const arma::vec& W) {
  double rss_local = arma::sum(W % arma::pow(y - yhat, 2));
  double tss_local = arma::sum(W % arma::pow(y - arma::mean(y), 2));
  return 1 - (rss_local / tss_local);
}

// [[Rcpp::export]]
Rcpp::NumericVector SelectSortedBW(const arma::mat& dist_mat,
                                   double start_idx, double end_idx) {
    int n_rows = dist_mat.n_rows;
    int n_cols = dist_mat.n_cols;
    int total_elements = n_rows * n_cols;
    start_idx = total_elements * start_idx;
    end_idx = total_elements * end_idx;

    arma::vec selected(end_idx - start_idx);
    arma::vec flat_dist_mat = arma::vectorise(dist_mat);
    selected = flat_dist_mat.subvec(start_idx, end_idx - 1);
    selected = arma::sort(selected);
    selected = arma::unique(selected);
    return ArmaVec4RcppNumericVector(selected);
}

arma::mat MatRowStandardize(const arma::mat& mat) {
  arma::mat result = mat;
  int n = result.n_rows;
  for (int i = 0; i < n; ++i) {
    double row_sum = arma::sum(result.row(i));
    if (row_sum != 0) {
      result.row(i) /= row_sum;
    }
  }
  return result;
}

arma::mat MatGlobalStandardize(const arma::mat& mat) {
  arma::mat result = mat;
  double mat_sum = arma::accu(result);
  if (mat_sum != 0) {
    result /= mat_sum;
  }
  return result;
}

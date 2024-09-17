#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Gaussian Kernel
double gaussian_kernel(double dist, double bw) {
  return exp(-0.5 * pow(dist / bw, 2));
}

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

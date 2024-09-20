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

double min_nona(NumericVector x) {
  NumericVector x1 = x[!is_na(x)];
  double xres = min(x1);
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

double CosineSimilarity(NumericVector A, NumericVector B) {
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

IntegerVector rcpp_which(LogicalVector x) {
  IntegerVector v = seq(0,x.size()-1);
  return v[x];
}

IntegerVector rcpp_seq(int start, int end) {
  int size = end - start + 1;
  IntegerVector result(size);
  for(int i = 0; i < size; ++i) {
    result[i] = start + i;
  }
  return result;
}

IntegerVector extract_window(IntegerMatrix mat, int window_size) {
  int rows = mat.nrow();
  int cols = mat.ncol();
  int half_window = window_size / 2;
  int result_size = rows * cols * window_size * window_size;
  IntegerVector result(result_size, R_NaInt);

  int index = 0;
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      for (int wr = -half_window; wr <= half_window; ++wr) {
        for (int wc = -half_window; wc <= half_window; ++wc) {
          int new_r = r + wr;
          int new_c = c + wc;
          if (new_r >= 0 && new_r < rows && new_c >= 0 && new_c < cols) {
            result[index] = mat(new_r, new_c);
          } else {
            result[index] = R_NaInt;
          }
          index++;
        }
      }
    }
  }
  return result;
}

NumericVector rcpp_log2(NumericVector vec) {
  NumericVector res(vec.size());
  for (int i = 0; i < vec.size(); ++i) {
    res[i] = std::log2(vec[i]);
  }
  return res;
}

NumericVector multiply_vector(NumericVector numvec1, NumericVector numvec2) {
  int n = numvec1.size();
  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = numvec1[i] * numvec2[i];
  }
  return result;
}

NumericVector NormalizeVector(NumericVector x, double a, double b) {
  double min_x = min(x);
  double max_x = max(x);
  if (min_x == max_x) {
    return NumericVector(x.size(), (a + b) / 2.0);
  }
  NumericVector res = a + (x - min_x) * (b - a) / (max_x - min_x);
  return res;
}

NumericMatrix NormalizeMatRow(NumericMatrix mat, double a, double b) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  for (int i = 0; i < nrow; ++i) {
    double row_min = mat(i, 0);
    double row_max = mat(i, 0);

    for (int j = 1; j < ncol; ++j) {
      if (mat(i, j) < row_min) row_min = mat(i, j);
      if (mat(i, j) > row_max) row_max = mat(i, j);
    }

    if (row_max == row_min) {
      for (int j = 0; j < ncol; ++j) {
        mat(i, j) = a;
      }
    } else {
      for (int j = 0; j < ncol; ++j) {
        mat(i, j) = a + (b - a) * (mat(i, j) - row_min) / (row_max - row_min);
      }
    }
  }

  return mat;
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
double InforEntropy(NumericVector x) {
  double res = -1 * sum_nona(x * rcpp_log2(x));
  return res;
}

NumericMatrix MatRowStandardize(NumericMatrix mat) {
  for (int i = 0; i < mat.nrow(); ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < mat.ncol(); ++j) {
      row_sum += mat(i, j);
    }
    if (row_sum != 0) {
      for (int j = 0; j < mat.ncol(); ++j) {
        mat(i, j) /= row_sum;
      }
    }
  }
  return mat;
}

NumericMatrix MatGlobalStandardize(NumericMatrix mat) {
  double mat_sum = 0.0;
  for (int i = 0; i < mat.nrow(); ++i) {
    for (int j = 0; j < mat.ncol(); ++j) {
      mat_sum += mat(i, j);
    }
  }
  for (int i = 0; i < mat.nrow(); ++i) {
    for (int j = 0; j < mat.ncol(); ++j) {
      mat(i, j) /= mat_sum;
    }
  }
  return mat;
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

NumericVector GCS_Variance(NumericMatrix x,
                           NumericMatrix wt,
                           String method) {
  NumericVector out(x.nrow());
  NumericVector theta(x.ncol());
  for (int i = 0; i < x.ncol(); ++i){
    theta[i] = pow(sd_nona(x(_,i)),2);
  }
  for (int i = 0; i < x.nrow(); ++i){
    NumericVector zs(x.nrow());
    NumericMatrix Ej(x.nrow(),x.ncol());
    for (int j = 0; j < x.ncol(); ++j){
      double zi = x(i,j);
      NumericVector delta(x.nrow());
      for (int n = 0; n < x.nrow(); ++n){
        double zn = x(n,j);
        delta[n] = pow((zi-zn),2);
      }
      double deltav = pow(mean_nona(delta),1/2);
      for (int nej = 0; nej < x.nrow(); ++nej){
        double znej = x(nej,j);
        Ej(nej,j) = exp(-pow((zi-znej),2)/(2*pow(theta[j]/deltav,2)));
      }
    }
    for (int ni = 0; ni < x.nrow(); ++ni){
      NumericVector zsn = Ej(ni,_);
      zs[ni] = min_nona(zsn);
    }
    if (method == "spvar"){
      out[i] = spatial_variance(zs,wt);
    } else {
      out[i] = InforEntropy(zs);
    }
  }
  return out;
}

NumericMatrix subset_matrix(NumericMatrix matrix,
                            IntegerVector rows,
                            IntegerVector cols) {
  NumericMatrix result(rows.size(), cols.size());
  for (int i = 0; i < rows.size(); ++i) {
    for (int j = 0; j < cols.size(); ++j) {
      result(i, j) = matrix(rows[i], cols[j]);
    }
  }
  return result;
}

NumericVector SSH_Variance(NumericVector x,
                           NumericMatrix wt,
                           String method) {
  NumericVector out(x.size());
  for (int i = 0; i < x.size(); ++i){
    IntegerVector wti_indice = rcpp_which(wt(i,_) != 0);
    NumericMatrix wti = subset_matrix(wt,wti_indice,wti_indice);
    NumericVector xi = x[wti_indice];
    if (method == "spvar"){
      out[i] = spatial_variance(xi,wt);
    } else {
      out[i] = InforEntropy(xi);
    }
  }
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

List RasterQueenNeighbors(int rows, int cols) {
  List neighbors(rows * cols);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      IntegerVector adj;
      for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
          if (di == 0 && dj == 0) continue;
          int ni = i + di;
          int nj = j + dj;
          if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
            adj.push_back(ni * cols + nj);
          }
        }
      }
      neighbors[i * cols + j] = adj;
    }
  }
  return neighbors;
}

List RasterGeoCNeighbors(int rows, int cols) {
  int center_row = rows / 2;
  int center_col = cols / 2;
  int center_idx = center_row * cols + center_col;

  List neighbors = RasterQueenNeighbors(rows, cols);
  IntegerVector center_neighbors = neighbors[center_idx];

  List result(rows*cols);
  for (int i = 0; i < rows; ++i) {
    IntegerVector row_intersection(cols);
    for (int j = 0; j < cols; ++j) {
      IntegerVector current_neighbors = neighbors[i * cols + j];
      IntegerVector intersection = intersect(center_neighbors, current_neighbors);
      result[i * rows + j] = intersection;
    }
  }
  return remove_index(result,center_idx);
}

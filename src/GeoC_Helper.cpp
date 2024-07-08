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

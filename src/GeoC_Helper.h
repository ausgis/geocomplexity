#ifndef GEO_C_HELPER_H
#define GEO_C_HELPER_H

#include <Rcpp.h>

double sum_nona(Rcpp::NumericVector x);
double mean_nona(Rcpp::NumericVector x);
double sd_nona(Rcpp::NumericVector x);
double matrix_sum(Rcpp::NumericMatrix mat);
double CosineSimilarity(Rcpp::NumericVector A, Rcpp::NumericVector B);
double spatial_variance(Rcpp::NumericVector x, Rcpp::NumericMatrix wt);
Rcpp::IntegerVector rcpp_which(Rcpp::LogicalVector x);
Rcpp::IntegerVector rcpp_seq(int start, int end);
Rcpp::IntegerVector extract_window(Rcpp::IntegerMatrix mat, int window_size);
Rcpp::NumericVector multiply_vector(Rcpp::NumericVector numvec1, Rcpp::NumericVector numvec2);
Rcpp::NumericVector GCS_Variance(Rcpp::NumericMatrix x, Rcpp::NumericMatrix wt);
Rcpp::List remove_index(Rcpp::List lst, int idx);
Rcpp::List RasterQueenNeighbors(int rows, int cols);
Rcpp::List RasterGeoCNeighbors(int rows, int cols);

#endif // GEO_C_HELPER_H

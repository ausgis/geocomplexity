#ifndef GEO_C_HELPER_H
#define GEO_C_HELPER_H

#include <Rcpp.h>

double sum_nona(Rcpp::NumericVector x);
double mean_nona(Rcpp::NumericVector x);
double sd_nona(Rcpp::NumericVector x);
double matrix_sum(Rcpp::NumericMatrix mat);
double cosine_similarity(double A, double B);
Rcpp::IntegerVector rcpp_which(Rcpp::LogicalVector x);
Rcpp::NumericVector multiply_vector(Rcpp::NumericVector numvec1, Rcpp::NumericVector numvec2);
double spatial_variance(Rcpp::NumericVector x, Rcpp::NumericMatrix wt);
Rcpp::List remove_index(Rcpp::List lst, int idx);

#endif // GEO_C_HELPER_H

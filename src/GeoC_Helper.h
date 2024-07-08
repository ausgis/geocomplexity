#ifndef GEO_C_HELPER_H
#define GEO_C_HELPER_H

#include <Rcpp.h>

double sum_nona(Rcpp::NumericVector x);
double mean_nona(Rcpp::NumericVector x);
Rcpp::IntegerVector rcpp_which(Rcpp::LogicalVector x);
Rcpp::NumericVector multiply_vector(Rcpp::NumericVector numvec1, Rcpp::NumericVector numvec2);
Rcpp::List remove_index(Rcpp::List lst, int idx);

#endif // GEO_C_HELPER_H

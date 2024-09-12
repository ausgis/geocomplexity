#ifndef GEOCGWR_HELPER_H
#define GEOCGWR_HELPER_H

#include <Rcpp.h>

double gaussian_kernel(double dist, double bw);
double exponential_kernel(double dist, double bw);
double bisquare_kernel(double dist, double bw);
double triangular_kernel(double dist, double bw);
double boxcar_kernel(double dist, double bw);

#endif // GEOCGWR_HELPER_H

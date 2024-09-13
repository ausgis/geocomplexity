#ifndef GEOCGWR_HELPER_H
#define GEOCGWR_HELPER_H

#include <RcppArmadillo.h>

double gaussian_kernel(double dist, double bw);
double exponential_kernel(double dist, double bw);
double bisquare_kernel(double dist, double bw);
double triangular_kernel(double dist, double bw);
double boxcar_kernel(double dist, double bw);
arma::vec Normalize4Interval(const arma::vec& v, double a, double b);
arma::mat DiagMatrix(int n);
arma::vec ArmaSeq(double start, double end, double step);

#endif // GEOCGWR_HELPER_H

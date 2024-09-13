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
arma::vec Double4Vec(double x);
arma::vec GenAdaptiveKNNBW(const arma::mat& D, double k);
double MaxInMatrix(const arma::mat& mat);
double MinInMatrix(const arma::mat& mat);

#endif // GEOCGWR_HELPER_H

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
arma::vec ArmaSeq(double from, double to, double by = 1, int length_out = -1);
arma::vec Double4Vec(double x);
arma::vec GenAdaptiveKNNBW(const arma::mat& D, double k);
double MaxInMatrix(const arma::mat& mat);
double MinInMatrix(const arma::mat& mat);
arma::mat EucdistM(const arma::mat& X);
double RowDiffAbsMean(const arma::mat& X, int i, int j);
double MeanWeight(const arma::vec& x, const arma::vec& w);

#endif // GEOCGWR_HELPER_H

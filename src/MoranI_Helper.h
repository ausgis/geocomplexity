#ifndef MORAN_I_HELPER_H
#define MORAN_I_HELPER_H

#include <Rcpp.h>
#include <RcppArmadillo.h>

std::string star(double p);
double pfunc(double z, std::string alternative);

#endif // MORAN_I_HELPER_H

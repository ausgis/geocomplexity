#include <Rcpp.h>
using namespace Rcpp;

// Gaussian Kernel
double gaussian_kernel(double dist, double bw) {
  return exp(-0.5 * pow(dist / bw, 2));
}

// Exponential Kernel
double exponential_kernel(double dist, double bw) {
  return exp(-dist / bw);
}

// Bisquare Kernel
double bisquare_kernel(double dist, double bw) {
  if (dist < bw) {
    return pow(1 - pow(dist / bw, 2), 2);
  } else {
    return 0.0;
  }
}

// Triangular Kernel
double triangular_kernel(double dist, double bw) {
  if (dist < bw) {
    return pow(1 - pow(dist / bw, 3), 3);
  } else {
    return 0.0;
  }
}

// Boxcar Kernel
double boxcar_kernel(double dist, double bw) {
  if (dist < bw) {
    return 1.0;
  } else {
    return 0.0;
  }
}

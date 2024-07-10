#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to calculate the star based on p-value
std::string star(double p) {
  if (p < 0.001) return "***";
  if (p < 0.01) return "**";
  if (p < 0.05) return "*";
  if (p < 0.1) return ".";
  return " ";
}

// Function to calculate p-value based on z and alternative hypothesis
double pfunc(double z, std::string alternative) {
  if (alternative == "greater") {
    return R::pnorm(z, 0.0, 1.0, false, false);
  } else if (alternative == "less") {
    return R::pnorm(z, 0.0, 1.0, true, false);
  } else { // "two.sided"
    return 2 * R::pnorm(std::abs(z), 0.0, 1.0, false, false);
  }
}

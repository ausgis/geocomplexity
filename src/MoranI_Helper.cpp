#include <Rcpp.h>
#include <iomanip>
#include <sstream>
#include <string>

using namespace Rcpp;

// Function to calculate the star based on p-value
std::string star(double p) {
  if (p < 0.001) return "***";
  if (p < 0.01) return "**";
  if (p < 0.05) return "*";
  if (p < 0.1) return ".";
  return " ";
}

// Helper function to concatenate a number and a significance star
std::string concat_num_with_star(double num, double p) {
  std::ostringstream oss;
  oss << num << star(p);
  return oss.str();
}

// Function to print the global spatial autocorrelation test results
// [[Rcpp::export]]
DataFrame print_global_moranI(DataFrame df) {
  CharacterVector variable = df["Variable"];
  NumericVector moran_i = df["MoranI"];
  NumericVector ei = df["EI"];
  NumericVector vari = df["VarI"];
  NumericVector zi = df["ZI"];
  NumericVector pi = df["PI"];
  CharacterVector stars(variable.size());

  for (int i = 0; i < variable.size(); ++i) {
    stars[i] = concat_num_with_star(moran_i[i],pi[i]);
  }

  DataFrame out = Rcpp::DataFrame::create(Rcpp::Named("Variable") = variable,
                                          Rcpp::Named("MoranI") = stars,
                                          Rcpp::Named("EI") = ei,
                                          Rcpp::Named("VarI") = vari,
                                          Rcpp::Named("zI") = zi,
                                          Rcpp::Named("pI") = pi);

  return out;
}

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

// Helper function to format a value to string with a fixed width
template <typename T>
std::string format_value(T value, int width, bool scientific = false) {
  std::ostringstream ss;
  if (scientific) {
    ss << std::scientific << std::setprecision(3);
  } else {
    ss << std::fixed << std::setprecision(4);
  }
  ss << std::left << std::setw(width) << value;
  return ss.str();
}

// Function to print the global spatial autocorrelation test results
// [[Rcpp::export]]
void print_global_moranI(DataFrame df) {
  std::vector<std::string> variable = df["Variable"];
  NumericVector moran_i = df["I"];
  NumericVector ei = df["EI"];
  NumericVector vari = df["VarI"];
  NumericVector zi = df["zI"];
  NumericVector pi = df["pI"];

  int width_var = 10;
  int width_val = 12;

  Rcout << " ***               global spatial autocorrelation test                  \n";
  Rcout << " -----------------------------------------------------------------------\n";
  Rcout << "   Variable      Moran I        EI         VarI       zI          pI    \n";
  Rcout << " ------------ ------------- ----------- ---------- --------- -----------\n";

  for (size_t i = 0; i < variable.size(); ++i) {
    std::string row = "   " + format_value(variable[i], width_val + 1) +
      star(pi[i]) + format_value(moran_i[i], width_var + 1) +
      format_value(ei[i], width_val) +
      format_value(vari[i], width_var) +
      format_value(zi[i], width_var) +
      format_value(pi[i], width_var, true);
    Rcout << row << "\n";
  }

  Rcout << " -----------------------------------------------------------------------\n";
}

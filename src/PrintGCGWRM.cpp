#include <Rcpp.h>
using namespace Rcpp;

void PrintTableMd(CharacterMatrix table, std::string indent = "  ") {
  for (int i = 0; i < table.nrow(); i++) {
    for (int j = 0; j < table.ncol(); j++) {
      Rcpp::Rcout << indent << as<std::string>(table(i, j));
      if (j < table.ncol() - 1) Rcpp::Rcout << " ";
    }
    Rcpp::Rcout << std::endl;
  }
}

Rcpp::CharacterMatrix Matrix2Char(NumericMatrix mat, std::string fmt = "%.3f") {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  CharacterMatrix result(nrow, ncol);

  char buffer[100];
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      std::sprintf(buffer, fmt.c_str(), mat(i, j));
      result(i, j) = std::string(buffer);
    }
  }
  return result;
}


// [[Rcpp::export]]
void PrintGCGWRM(Rcpp::List x, Rcpp::DataFrame coefs) {
  // Basic Information
  Rcpp::Rcout << "Geographical Complexity-Geographically Weighted Regression Model" << std::endl;
  Rcpp::Rcout << "================================================================" << std::endl;
  Rcpp::Rcout << std::endl;

  // Diagnostic Information
  Rcpp::Rcout << "Diagnostic Information" << std::endl;
  Rcpp::Rcout << "----------------------" << std::endl;
  List diagnostic = x["diagnostic"];
  Rcpp::Rcout << "  RSS: " << as<double>(diagnostic["RSS"]) << std::endl;
  Rcpp::Rcout << "  ENP: " << as<double>(diagnostic["ENP"]) << std::endl;
  Rcpp::Rcout << "  EDF: " << as<double>(diagnostic["EDF"]) << std::endl;
  Rcpp::Rcout << "   R2: " << as<double>(diagnostic["R2"]) << std::endl;
  Rcpp::Rcout << "R2adj: " << as<double>(diagnostic["R2_Adj"]) << std::endl;
  Rcpp::Rcout << "  AIC: " << as<double>(diagnostic["AIC"]) << std::endl;
  Rcpp::Rcout << " AICc: " << as<double>(diagnostic["AICc"]) << std::endl;
  Rcpp::Rcout << " RMSE: " << as<double>(diagnostic["RMSE"]) << std::endl;
  Rcpp::Rcout << std::endl;

  // Summary of Coefficient Estimates
  Rcpp::Rcout << "Summary of Coefficient Estimates" << std::endl;
  Rcpp::Rcout << "--------------------------------" << std::endl;
}

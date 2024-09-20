#include <Rcpp.h>
using namespace Rcpp;

std::string FormatBW(int knn, const std::string& criterion) {
  char buffer[100];
  std::sprintf(buffer, "%d (Nearest Neighbours) (Optimized according to %s)", knn, criterion.c_str());
  return std::string(buffer);
}

// [[Rcpp::export]]
void PrintGCGWRM(Rcpp::List x) {
  // Basic Information
  Rcpp::Rcout << "Geographical Complexity-Geographically Weighted Regression Model" << std::endl;
  Rcpp::Rcout << "================================================================" << std::endl;
  Rcpp::List args = x["args"];
  bool adaptive = x["adaptive"];
  int knn = args["knn"];
  double alpha = args["alpha"];
  std::string kernel = args["kernel"];
  std::string criterion = args["criterion"];
  std::string bw = as<std::string>(args["bw"]);
  if (adaptive) {
    bw = FormatBW(knn,criterion);
  }
  Rcpp::Rcout << "     Kernel:  " <<  kernel << std::endl;
  Rcpp::Rcout << "  Bandwidth:  " <<  bw     << std::endl;
  Rcpp::Rcout << "      Alpha:  " <<  alpha  << std::endl;
  Rcpp::Rcout << std::endl;

  // Diagnostic Information
  Rcpp::Rcout << "Diagnostic Information" << std::endl;
  Rcpp::Rcout << "----------------------" << std::endl;
  Rcpp::List diagnostic = x["diagnostic"];
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

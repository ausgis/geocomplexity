#include <Rcpp.h>
#include <iomanip>
using namespace Rcpp;

std::string FormatBW(int knn, const std::string& criterion) {
  char buffer[100];
  std::sprintf(buffer, "%d (Nearest Neighbours) (Optimized according to %s)", knn, criterion.c_str());
  return std::string(buffer);
}

void PrintCoefMat(Rcpp::NumericMatrix mat, Rcpp::CharacterVector rownames) {
  // Set up the width for columns, similar to the R print table format
  int colWidth = 10;
  Rcpp::CharacterVector colnames = Rcpp::CharacterVector::create("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.");
  Rcpp::Rcout << std::left << std::setw(colWidth) << "Coefficient";

  // Print column headers
  for (int j = 0; j < colnames.size(); ++j) {
    Rcpp::Rcout << std::setw(colWidth) << colnames[j];
  }
  Rcpp::Rcout << std::endl;

  // Print the matrix with row names
  for (int i = 0; i < mat.nrow(); ++i) {
    Rcpp::Rcout << std::left << std::setw(colWidth) << rownames[i];
    for (int j = 0; j < mat.ncol(); ++j) {
      Rcpp::Rcout << std::setw(colWidth) << std::setprecision(3) << std::fixed << mat(i, j);
    }
    Rcpp::Rcout << std::endl;
  }
  Rcpp::Rcout << std::endl;
}

// [[Rcpp::export]]
void PrintGCGWRM(Rcpp::List x,
                 Rcpp::NumericMatrix coefmat,
                 Rcpp::CharacterVector coefname) {
  // Basic Information
  Rcpp::Rcout << "Geographical Complexity-Geographically Weighted Regression Model" << std::endl;
  Rcpp::Rcout << "================================================================" << std::endl;
  Rcpp::List args = x["args"];
  bool adaptive = args["adaptive"];
  int knn = args["knn"];
  double alpha = args["alpha"];
  std::string kernel = args["kernel"];
  std::string criterion = args["criterion"];
  std::string bwstr = FormatBW(knn,criterion);
  double bw1 = args["bw"];
  std::string bw = std::to_string(bw1);

  Rcpp::Rcout << "     Kernel:  " <<  kernel << std::endl;
  Rcpp::Rcout << "  Bandwidth:  " <<  (adaptive ? bwstr : bw)  << std::endl;
  Rcpp::Rcout << "      Alpha:  " <<  alpha  << std::endl;
  Rcpp::Rcout << std::endl;

  // Summary of Coefficient Estimates
  Rcpp::Rcout << "Summary of Coefficient Estimates" << std::endl;
  Rcpp::Rcout << "--------------------------------" << std::endl;
  PrintCoefMat(coefmat,coefname);

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
}

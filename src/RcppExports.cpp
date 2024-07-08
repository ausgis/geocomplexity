// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// spatial_variance
double spatial_variance(NumericVector x, NumericMatrix wt);
RcppExport SEXP _geocomplexity_spatial_variance(SEXP xSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_variance(x, wt));
    return rcpp_result_gen;
END_RCPP
}
// RasterQueenNeighbors
List RasterQueenNeighbors(int rows, int cols);
RcppExport SEXP _geocomplexity_RasterQueenNeighbors(SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterQueenNeighbors(rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCNeighbors
List RasterGeoCNeighbors(int rows, int cols);
RcppExport SEXP _geocomplexity_RasterGeoCNeighbors(SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCNeighbors(rows, cols));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCDependenceOne
double RasterGeoCDependenceOne(NumericVector x, size_t ni, size_t nw);
RcppExport SEXP _geocomplexity_RasterGeoCDependenceOne(SEXP xSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCDependenceOne(x, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCDependence
NumericVector RasterGeoCDependence(NumericVector x, size_t ni, size_t nw);
RcppExport SEXP _geocomplexity_RasterGeoCDependence(SEXP xSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCDependence(x, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// VectorGeoCDependence
NumericVector VectorGeoCDependence(NumericVector x, NumericMatrix wt);
RcppExport SEXP _geocomplexity_VectorGeoCDependence(SEXP xSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorGeoCDependence(x, wt));
    return rcpp_result_gen;
END_RCPP
}
// VectorGeoCSimilarity
NumericVector VectorGeoCSimilarity(NumericVector x, NumericMatrix wt);
RcppExport SEXP _geocomplexity_VectorGeoCSimilarity(SEXP xSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorGeoCSimilarity(x, wt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geocomplexity_spatial_variance", (DL_FUNC) &_geocomplexity_spatial_variance, 2},
    {"_geocomplexity_RasterQueenNeighbors", (DL_FUNC) &_geocomplexity_RasterQueenNeighbors, 2},
    {"_geocomplexity_RasterGeoCNeighbors", (DL_FUNC) &_geocomplexity_RasterGeoCNeighbors, 2},
    {"_geocomplexity_RasterGeoCDependenceOne", (DL_FUNC) &_geocomplexity_RasterGeoCDependenceOne, 3},
    {"_geocomplexity_RasterGeoCDependence", (DL_FUNC) &_geocomplexity_RasterGeoCDependence, 3},
    {"_geocomplexity_VectorGeoCDependence", (DL_FUNC) &_geocomplexity_VectorGeoCDependence, 2},
    {"_geocomplexity_VectorGeoCSimilarity", (DL_FUNC) &_geocomplexity_VectorGeoCSimilarity, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_geocomplexity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

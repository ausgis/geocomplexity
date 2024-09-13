// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BasicGWRFit
Rcpp::List BasicGWRFit(arma::vec y, arma::mat X, arma::mat Cdist, double bw, double knn, bool adaptive, std::string kernel);
RcppExport SEXP _geocomplexity_BasicGWRFit(SEXP ySEXP, SEXP XSEXP, SEXP CdistSEXP, SEXP bwSEXP, SEXP knnSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cdist(CdistSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< double >::type knn(knnSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(BasicGWRFit(y, X, Cdist, bw, knn, adaptive, kernel));
    return rcpp_result_gen;
END_RCPP
}
// BasicGWRSel
Rcpp::List BasicGWRSel(arma::vec bandwidth, arma::vec knns, arma::vec y, arma::mat X, arma::mat Cdist, bool adaptive, std::string criterion, std::string kernel);
RcppExport SEXP _geocomplexity_BasicGWRSel(SEXP bandwidthSEXP, SEXP knnsSEXP, SEXP ySEXP, SEXP XSEXP, SEXP CdistSEXP, SEXP adaptiveSEXP, SEXP criterionSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type knns(knnsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cdist(CdistSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< std::string >::type criterion(criterionSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(BasicGWRSel(bandwidth, knns, y, X, Cdist, adaptive, criterion, kernel));
    return rcpp_result_gen;
END_RCPP
}
// BasicGWR
Rcpp::List BasicGWR(arma::vec y, arma::mat X, arma::mat Cdist, SEXP bw, bool adaptive, std::string kernel);
RcppExport SEXP _geocomplexity_BasicGWR(SEXP ySEXP, SEXP XSEXP, SEXP CdistSEXP, SEXP bwSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cdist(CdistSEXP);
    Rcpp::traits::input_parameter< SEXP >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(BasicGWR(y, X, Cdist, bw, adaptive, kernel));
    return rcpp_result_gen;
END_RCPP
}
// GeoCGWRFit
Rcpp::List GeoCGWRFit(arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist, double bwc, double bwg, double knn, bool adaptive, double alpha, std::string kernel);
RcppExport SEXP _geocomplexity_GeoCGWRFit(SEXP ySEXP, SEXP XSEXP, SEXP gcsSEXP, SEXP GdistSEXP, SEXP bwcSEXP, SEXP bwgSEXP, SEXP knnSEXP, SEXP adaptiveSEXP, SEXP alphaSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gcs(gcsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gdist(GdistSEXP);
    Rcpp::traits::input_parameter< double >::type bwc(bwcSEXP);
    Rcpp::traits::input_parameter< double >::type bwg(bwgSEXP);
    Rcpp::traits::input_parameter< double >::type knn(knnSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(GeoCGWRFit(y, X, gcs, Gdist, bwc, bwg, knn, adaptive, alpha, kernel));
    return rcpp_result_gen;
END_RCPP
}
// GeoCGWRSel
Rcpp::List GeoCGWRSel(arma::vec bwcs, arma::vec bwgs, arma::vec knns, arma::vec alpha, arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist, bool adaptive, std::string criterion, std::string kernel);
RcppExport SEXP _geocomplexity_GeoCGWRSel(SEXP bwcsSEXP, SEXP bwgsSEXP, SEXP knnsSEXP, SEXP alphaSEXP, SEXP ySEXP, SEXP XSEXP, SEXP gcsSEXP, SEXP GdistSEXP, SEXP adaptiveSEXP, SEXP criterionSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bwcs(bwcsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bwgs(bwgsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type knns(knnsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gcs(gcsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gdist(GdistSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< std::string >::type criterion(criterionSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(GeoCGWRSel(bwcs, bwgs, knns, alpha, y, X, gcs, Gdist, adaptive, criterion, kernel));
    return rcpp_result_gen;
END_RCPP
}
// GeoCGWR
Rcpp::List GeoCGWR(arma::vec y, arma::mat X, arma::mat gcs, arma::mat Gdist, SEXP bwc, SEXP bwg, arma::vec alpha, bool adaptive, std::string kernel);
RcppExport SEXP _geocomplexity_GeoCGWR(SEXP ySEXP, SEXP XSEXP, SEXP gcsSEXP, SEXP GdistSEXP, SEXP bwcSEXP, SEXP bwgSEXP, SEXP alphaSEXP, SEXP adaptiveSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gcs(gcsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gdist(GdistSEXP);
    Rcpp::traits::input_parameter< SEXP >::type bwc(bwcSEXP);
    Rcpp::traits::input_parameter< SEXP >::type bwg(bwgSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(GeoCGWR(y, X, gcs, Gdist, bwc, bwg, alpha, adaptive, kernel));
    return rcpp_result_gen;
END_RCPP
}
// InforEntropy
double InforEntropy(NumericVector x);
RcppExport SEXP _geocomplexity_InforEntropy(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(InforEntropy(x));
    return rcpp_result_gen;
END_RCPP
}
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
// GeoCSWT
NumericMatrix GeoCSWT(NumericVector x, NumericMatrix wt, String style);
RcppExport SEXP _geocomplexity_GeoCSWT(SEXP xSEXP, SEXP wtSEXP, SEXP styleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< String >::type style(styleSEXP);
    rcpp_result_gen = Rcpp::wrap(GeoCSWT(x, wt, style));
    return rcpp_result_gen;
END_RCPP
}
// MI_vec
Rcpp::DataFrame MI_vec(arma::mat x, arma::mat W, std::string alternative, bool symmetrize);
RcppExport SEXP _geocomplexity_MI_vec(SEXP xSEXP, SEXP WSEXP, SEXP alternativeSEXP, SEXP symmetrizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type symmetrize(symmetrizeSEXP);
    rcpp_result_gen = Rcpp::wrap(MI_vec(x, W, alternative, symmetrize));
    return rcpp_result_gen;
END_RCPP
}
// print_global_moranI
DataFrame print_global_moranI(DataFrame df);
RcppExport SEXP _geocomplexity_print_global_moranI(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(print_global_moranI(df));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCMoranOne
double RasterGeoCMoranOne(NumericVector x, size_t ni, size_t nw);
RcppExport SEXP _geocomplexity_RasterGeoCMoranOne(SEXP xSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCMoranOne(x, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCMoran
NumericVector RasterGeoCMoran(NumericVector x, size_t ni, size_t nw);
RcppExport SEXP _geocomplexity_RasterGeoCMoran(SEXP xSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCMoran(x, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCSSH
NumericVector RasterGeoCSSH(NumericVector x, IntegerMatrix iw, int w, String method);
RcppExport SEXP _geocomplexity_RasterGeoCSSH(SEXP xSEXP, SEXP iwSEXP, SEXP wSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type iw(iwSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCSSH(x, iw, w, method));
    return rcpp_result_gen;
END_RCPP
}
// RasterGeoCSimilarity
NumericVector RasterGeoCSimilarity(NumericMatrix x, IntegerMatrix iw, int w, int similarity, String method);
RcppExport SEXP _geocomplexity_RasterGeoCSimilarity(SEXP xSEXP, SEXP iwSEXP, SEXP wSEXP, SEXP similaritySEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type iw(iwSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(RasterGeoCSimilarity(x, iw, w, similarity, method));
    return rcpp_result_gen;
END_RCPP
}
// VectorGeoCMoran
NumericVector VectorGeoCMoran(NumericVector x, NumericMatrix wt);
RcppExport SEXP _geocomplexity_VectorGeoCMoran(SEXP xSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorGeoCMoran(x, wt));
    return rcpp_result_gen;
END_RCPP
}
// VectorGeoCSSH
NumericVector VectorGeoCSSH(NumericVector xobs, NumericMatrix wt, String method);
RcppExport SEXP _geocomplexity_VectorGeoCSSH(SEXP xobsSEXP, SEXP wtSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorGeoCSSH(xobs, wt, method));
    return rcpp_result_gen;
END_RCPP
}
// VectorGeoCSimilarity
NumericVector VectorGeoCSimilarity(NumericMatrix xobs, NumericMatrix wt, int similarity, String method);
RcppExport SEXP _geocomplexity_VectorGeoCSimilarity(SEXP xobsSEXP, SEXP wtSEXP, SEXP similaritySEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorGeoCSimilarity(xobs, wt, similarity, method));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geocomplexity_BasicGWRFit", (DL_FUNC) &_geocomplexity_BasicGWRFit, 7},
    {"_geocomplexity_BasicGWRSel", (DL_FUNC) &_geocomplexity_BasicGWRSel, 8},
    {"_geocomplexity_BasicGWR", (DL_FUNC) &_geocomplexity_BasicGWR, 6},
    {"_geocomplexity_GeoCGWRFit", (DL_FUNC) &_geocomplexity_GeoCGWRFit, 10},
    {"_geocomplexity_GeoCGWRSel", (DL_FUNC) &_geocomplexity_GeoCGWRSel, 11},
    {"_geocomplexity_GeoCGWR", (DL_FUNC) &_geocomplexity_GeoCGWR, 9},
    {"_geocomplexity_InforEntropy", (DL_FUNC) &_geocomplexity_InforEntropy, 1},
    {"_geocomplexity_spatial_variance", (DL_FUNC) &_geocomplexity_spatial_variance, 2},
    {"_geocomplexity_GeoCSWT", (DL_FUNC) &_geocomplexity_GeoCSWT, 3},
    {"_geocomplexity_MI_vec", (DL_FUNC) &_geocomplexity_MI_vec, 4},
    {"_geocomplexity_print_global_moranI", (DL_FUNC) &_geocomplexity_print_global_moranI, 1},
    {"_geocomplexity_RasterGeoCMoranOne", (DL_FUNC) &_geocomplexity_RasterGeoCMoranOne, 3},
    {"_geocomplexity_RasterGeoCMoran", (DL_FUNC) &_geocomplexity_RasterGeoCMoran, 3},
    {"_geocomplexity_RasterGeoCSSH", (DL_FUNC) &_geocomplexity_RasterGeoCSSH, 4},
    {"_geocomplexity_RasterGeoCSimilarity", (DL_FUNC) &_geocomplexity_RasterGeoCSimilarity, 5},
    {"_geocomplexity_VectorGeoCMoran", (DL_FUNC) &_geocomplexity_VectorGeoCMoran, 2},
    {"_geocomplexity_VectorGeoCSSH", (DL_FUNC) &_geocomplexity_VectorGeoCSSH, 3},
    {"_geocomplexity_VectorGeoCSimilarity", (DL_FUNC) &_geocomplexity_VectorGeoCSimilarity, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_geocomplexity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

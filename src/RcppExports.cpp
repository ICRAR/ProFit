// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// profitMakeSersic
NumericMatrix profitMakeSersic(const IntegerMatrix& CALCREGION, const double XCEN, const double YCEN, const double MAG, const double RE, const double NSER, const double ANG, const double AXRAT, const double BOX, const double MAGZERO, const bool ROUGH, const NumericVector& XLIM, const NumericVector& YLIM, const IntegerVector& DIM, const int UPSCALE, const int MAXDEPTH, const double RESWITCH, const double ACC, const bool DOCALCREGION);
RcppExport SEXP ProFit_profitMakeSersic(SEXP CALCREGIONSEXP, SEXP XCENSEXP, SEXP YCENSEXP, SEXP MAGSEXP, SEXP RESEXP, SEXP NSERSEXP, SEXP ANGSEXP, SEXP AXRATSEXP, SEXP BOXSEXP, SEXP MAGZEROSEXP, SEXP ROUGHSEXP, SEXP XLIMSEXP, SEXP YLIMSEXP, SEXP DIMSEXP, SEXP UPSCALESEXP, SEXP MAXDEPTHSEXP, SEXP RESWITCHSEXP, SEXP ACCSEXP, SEXP DOCALCREGIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type CALCREGION(CALCREGIONSEXP);
    Rcpp::traits::input_parameter< const double >::type XCEN(XCENSEXP);
    Rcpp::traits::input_parameter< const double >::type YCEN(YCENSEXP);
    Rcpp::traits::input_parameter< const double >::type MAG(MAGSEXP);
    Rcpp::traits::input_parameter< const double >::type RE(RESEXP);
    Rcpp::traits::input_parameter< const double >::type NSER(NSERSEXP);
    Rcpp::traits::input_parameter< const double >::type ANG(ANGSEXP);
    Rcpp::traits::input_parameter< const double >::type AXRAT(AXRATSEXP);
    Rcpp::traits::input_parameter< const double >::type BOX(BOXSEXP);
    Rcpp::traits::input_parameter< const double >::type MAGZERO(MAGZEROSEXP);
    Rcpp::traits::input_parameter< const bool >::type ROUGH(ROUGHSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type XLIM(XLIMSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type YLIM(YLIMSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type DIM(DIMSEXP);
    Rcpp::traits::input_parameter< const int >::type UPSCALE(UPSCALESEXP);
    Rcpp::traits::input_parameter< const int >::type MAXDEPTH(MAXDEPTHSEXP);
    Rcpp::traits::input_parameter< const double >::type RESWITCH(RESWITCHSEXP);
    Rcpp::traits::input_parameter< const double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< const bool >::type DOCALCREGION(DOCALCREGIONSEXP);
    __result = Rcpp::wrap(profitMakeSersic(CALCREGION, XCEN, YCEN, MAG, RE, NSER, ANG, AXRAT, BOX, MAGZERO, ROUGH, XLIM, YLIM, DIM, UPSCALE, MAXDEPTH, RESWITCH, ACC, DOCALCREGION));
    return __result;
END_RCPP
}
// profitBruteConv
NumericMatrix profitBruteConv(NumericMatrix image, NumericMatrix psf, const IntegerMatrix& CALCREGION, const bool DOCALCREGION);
RcppExport SEXP ProFit_profitBruteConv(SEXP imageSEXP, SEXP psfSEXP, SEXP CALCREGIONSEXP, SEXP DOCALCREGIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type image(imageSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type psf(psfSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type CALCREGION(CALCREGIONSEXP);
    Rcpp::traits::input_parameter< const bool >::type DOCALCREGION(DOCALCREGIONSEXP);
    __result = Rcpp::wrap(profitBruteConv(image, psf, CALCREGION, DOCALCREGION));
    return __result;
END_RCPP
}

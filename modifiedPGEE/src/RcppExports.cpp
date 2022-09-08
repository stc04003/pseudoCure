// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SHM_CppArma
Rcpp::List SHM_CppArma(arma::vec y, arma::vec mu, arma::mat X, arma::cube Rhat, int N, int nx, arma::vec nt, arma::vec index, arma::vec muEta);
RcppExport SEXP _modifiedPGEE_SHM_CppArma(SEXP ySEXP, SEXP muSEXP, SEXP XSEXP, SEXP RhatSEXP, SEXP NSEXP, SEXP nxSEXP, SEXP ntSEXP, SEXP indexSEXP, SEXP muEtaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Rhat(RhatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nt(ntSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muEta(muEtaSEXP);
    rcpp_result_gen = Rcpp::wrap(SHM_CppArma(y, mu, X, Rhat, N, nx, nt, index, muEta));
    return rcpp_result_gen;
END_RCPP
}

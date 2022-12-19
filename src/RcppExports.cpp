// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4BinCauchy_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4BinFullborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4BinNoborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4BinNormal_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4ContCauchy_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4ContFullborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4ContNoborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4ContNormal_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4T2ECauchy_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4T2EFullborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4T2ENoborrow_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4T2ENormal_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4BinCauchy_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BinCauchy_mod, 0},
    {"_rcpp_module_boot_stan_fit4BinFullborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BinFullborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4BinNoborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BinNoborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4BinNormal_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BinNormal_mod, 0},
    {"_rcpp_module_boot_stan_fit4ContCauchy_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4ContCauchy_mod, 0},
    {"_rcpp_module_boot_stan_fit4ContFullborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4ContFullborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4ContNoborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4ContNoborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4ContNormal_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4ContNormal_mod, 0},
    {"_rcpp_module_boot_stan_fit4T2ECauchy_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4T2ECauchy_mod, 0},
    {"_rcpp_module_boot_stan_fit4T2EFullborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4T2EFullborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4T2ENoborrow_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4T2ENoborrow_mod, 0},
    {"_rcpp_module_boot_stan_fit4T2ENormal_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4T2ENormal_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_psBayesborrow(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
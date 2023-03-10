# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("BinCauchy", "BinFullborrow", "BinNoborrow", "BinNormal", "ContCauchy", "ContFullborrow", "ContNoborrow", "ContNormal", "T2ECauchy", "T2ECauchyC0", "T2EFullborrow", "T2EFullborrowC0", "T2ENoborrow", "T2ENoborrowC0", "T2ENormal", "T2ENormalC0")

# load each stan module
Rcpp::loadModule("stan_fit4BinCauchy_mod", what = TRUE)
Rcpp::loadModule("stan_fit4BinFullborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4BinNoborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4BinNormal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContCauchy_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContFullborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContNoborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContNormal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ECauchy_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ECauchyC0_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EFullborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EFullborrowC0_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ENoborrow_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ENoborrowC0_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ENormal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2ENormalC0_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})

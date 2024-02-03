
#' Simulation study of hybrid control design with Bayesian dynamic borrowing
#' incorporating propensity score matched external control: time-to-event outcome
#'
#' Simulation study is conducted to assess operating characteristics of hybrid
#' control design with Bayesian dynamic borrowing, where the concurrent control
#' is augmented by external control. The external controls are selected from
#' external control pool using a propensity score matching. Commensurate power
#' prior is used for Bayesian dynamic borrowing. The time-to-event outcome is
#' applicable.
#' @usage
#' psborrow.t2e(
#'   n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, n.EC, accrual,
#'   out.mevent.CT, out.mevent.CC, driftHR,
#'   cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
#'   psmatch.cov,
#'   method.psest="glm", method.pslink="logit",
#'   method.whomatch, method.matching, method.psorder, n.boot=100,
#'   analysis.cov, method.borrow,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025, nsim)
#' @param n.CT Number of patients in treatment group in the current trial.
#' @param n.CC Number of patients in concurrent control group in the current
#' trial.
#' @param nevent.C Number of events in treatment and concurrent control group
#' in the current trial.
#' @param n.ECp Number of patients in external control pool.
#' @param nevent.ECp Number of events in external control pool.
#' @param n.EC Number of patients in external control.
#' @param accrual Accrual rate, defined as the number of enrolled patients per
#' month.
#' @param out.mevent.CT True median time to event in treatment group in the
#' current trial.
#' @param out.mevent.CC True median time to event in concurrent control group
#' in the current trial.
#' @param driftHR Hazard ratio between concurrent and external control for
#' which the bias should be plotted.
#' @param cov.C List of covariate distributions for treatment and concurrent
#' control group in the current trial. Continuous and binary covariate are
#' applicable. The continuous covariate is assumed to follow a normal
#' distribution; for example, specified as
#' \code{list(dist="norm", mean=0, sd=1, lab="cov1")}. The binary covariate is
#' assumed to follow a binomial distribution; for example, specified as
#' \code{list(dist="binom", prob=0.4, lab="cov2")}. \code{lab} is the column
#' name of the covariate in the data frame generated.
#' @param cov.cor.C Matrix of correlation coefficients for each pair of
#' covariate for treatment and concurrent control group in the current trial,
#' specified as Gaussian copula parameter.
#' @param cov.EC List of covariate distributions for external control. The
#' continuous covariate is assumed to follow a normal distribution; for example,
#' specified as \code{list(dist="norm", mean=0, sd=1, lab="cov1")}. The binary
#' covariate is assumed to follow a binomial distribution; for example,
#' specified as \code{list(dist="binom", prob=0.4, lab="cov2")}. \code{lab} is
#' the column name of the covariate in the data frame generated, which must be
#' consistent with those used for \code{cov.C}.
#' @param cov.cor.EC Matrix of correlation coefficients for each pair of
#' covariate for external control, specified as Gaussian copula parameter.
#' @param cov.effect Vector of covariate effects on the outcome, specified as
#' hazard ratio per one unit increase in continuous covariate or as hazard
#' ratio between categories for binary covariate.
#' @param psmatch.cov Vector of names of covariates which are used for the
#' propensity score matching. The names of covariates must be included in
#' \code{lab} values specified in \code{cov.C}.
#' @param method.psest Method of estimating the propensity score. Allowable
#' options include, for example, \code{"glm"} for generalized linear model
#' (e.g., logistic regression); \code{"gam"} for generalized additive model;
#' \code{"gbm"} for generalized boosted model; \code{"lasso"} for lasso
#' regression; \code{"rpart"} for classification tree. The default value is
#' \code{method.psest="glm"}.
#' @param method.pslink Link function used in estimating the propensity score.
#' Allowable options depend on the specific \code{method.psest} value specified.
#' The default value is \code{method.pslink="logit"}, which, along with
#' \code{method.psest="glm"}, identifies the default method as logistic
#' regression.
#' @param method.whomatch Options of who to match. Allowable options include
#' \code{conc.contl} for matching concurrent control to external control pool;
#' \code{conc.treat} for matching treatment to external control pool;
#' \code{conc.all} for matching treatment plus concurrent control to external
#' control pool; \code{treat2contl} for matching treatment to concurrent control
#' plus external control pool.
#' @param method.matching Matching method. Allowable options include
#' \code{"optimal"} for optimal matching; \code{"nearest"} for nearest neighbor
#' matching without replacement; \code{"med.optimal"} for equally splitting
#' patients in the current trial and taking the median of each subset, followed
#' by 1:1 optimal matching; \code{"med.nearest"} for equally splitting
#' patients in the current trial and taking the median of each subset, followed
#' by 1:1 nearest neighbor matching without replacement; \code{"km.optimal"} for
#' k-means clustering of patients in the current trial, followed by 1:1 optimal
#' matching; \code{"km.nearest"} for k-means clustering of patients in the
#' current trial, followed by 1:1 nearest neighbor matching without replacement;
#' \code{"cm.optimal"} for fuzzy c-means clustering of patients in the current
#' trial, followed by 1:1 optimal matching; \code{"cm.nearest"} for fuzzy
#' c-means of patients in the current trial, followed by 1:1 nearest neighbor
#' matching without replacement; \code{"boot.optimal"} for bootstrap sampling
#' from patients in the current trial, followed by 1:1 optimal matching;
#' \code{"boot.nearest"} for bootstrap sampling from patient in the current
#' trial, followed by 1:1 nearest neighbor matching without replacement.
#' @param method.psorder Order that the matching takes place when a nearest
#' neighbor matching is used. Allowable options include \code{"largest"},
#' where matching takes place in descending order of propensity score;
#' \code{"smallest"}, where matching takes place in ascending order of
#' propensity score; \code{"random"}, where matching takes place in a random
#' order; \code{"data"}, where matching takes place based on the order of units
#' in the data. The matching order must be specified when using the nearest
#' neighbor matching.
#' @param n.boot Number of bootstrap sampling, which must be specified when
#' \code{method.matching="boot.optimal"} or
#' \code{method.matching="boot.nearest"}. The default value is \code{n.boot=100}.
#' @param analysis.cov Vector of names of covariates which are used for the
#' Bayesian analysis with commensurate prior. The names of covariates must be
#' included in \code{lab} values specified in \code{cov.C}.
#' @param method.borrow List of information borrowing method. \code{"noborrow"}
#' uses the concurrent data only. \code{"fullborrow"} uses the external control
#' data without discounting. \code{"cauchy"} uses the commensurate prior to
#' dynamically borrow the external control data, and the commensurability
#' parameter is assumed to follow a half-Cauchy distribution. \code{"normal"}
#' uses the commensurate prior to dynamically borrow the external control data,
#' and the commensurability parameter is assumed to follow a half-normal
#' distribution. \code{"cauchy"} and \code{"normal"} require to specify the
#' scale parameter \code{scale} of half-Cauchy and half-normal distribution
#' respectively.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=2}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=4000}.
#' @param warmup Number of warmup (burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.025}.
#' @param nsim Number of simulated trials.
#' @details The simulation study consists of three part: data generation
#' conducted by \code{trial.simulation.t2e} function, propensity score matching
#' conducted by \code{psmatch} function, and Bayesian analysis with commensurate
#' prior conducted by \code{commensurate.t2e} function. Users can specify
#' different sets of covariates for the propensity score matching and the
#' Bayesian analysis.
#' @return
#' The \code{psborrow.t2e} returns a list containing the following objects:
#' \item{reject}{Data frame containing results of Bayesian one-sided hypothesis
#' testing (whether or not the posterior probability that the log hazard ratio
#' is greater or less than 0 exceeds 1 minus significance level): \code{TRUE}
#' when significant, otherwise \code{FALSE}.}
#' \item{theta}{Data frame containing posterior mean, median, and sd of log
#' hazard ratio.}
#' \item{ov}{Data frame containing overlapping coefficient of propensity score
#' densities between treatment versus concurrent control plus external control}
#' \item{n.CT}{Number of patients in treatment group in the current trial.}
#' \item{n.CC}{Number of patients in concurrent control group in the current
#' trial.}
#' \item{n.ECp}{Number of patients in external control pool.}
#' \item{n.EC}{Number of patients in external control.}
#' \item{drift}{Hazard ratio between concurrent and external control.}
#' \item{true.theta}{True log hazard ratio}
#' \item{method.psest}{Method of estimating the propensity score.}
#' \item{method.pslink}{Link function used in estimating the propensity score.}
#' \item{method.whomatch}{Option of who to match.}
#' \item{method.matching}{Propensity score matching method.}
#' \item{method.psorder}{Order that the matching takes place when a nearest
#' neighbor matching is used.}
#' @examples
#' n.CT       <- 100
#' n.CC       <- 50
#' nevent.C   <- 100
#' n.ECp      <- 200
#' nevent.ECp <- 180
#' n.EC       <- 50
#' accrual    <- 16
#'
#' out.mevent.CT <- 6
#' out.mevent.CC <- 6
#' driftHR       <- 1
#'
#' cov.C <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'               list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.C <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'                list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.EC <- rbind(c(  1,0.1),
#'                     c(0.1,  1))
#'
#' cov.effect <- c(0.9,0.9)
#'
#' psmatch.cov <- c("cov1","cov2")
#'
#' method.whomatch <- "conc.treat"
#' method.matching <- "optimal"
#' method.psorder  <- NULL
#'
#' analysis.cov <- c("cov1")
#'
#' method.borrow <- list(list(prior="noborrow"),
#'                       list(prior="normal",scale=0.5))
#'
#' nsim <- 5
#'
#' psborrow.t2e(
#'   n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'   n.ECp=n.ECp, nevent.ECp=nevent.ECp, n.EC=n.EC, accrual=accrual,
#'   out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C,
#'   cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect,
#'   psmatch.cov=psmatch.cov, method.whomatch=method.whomatch,
#'   method.matching=method.matching, method.psorder=method.psorder,
#'   analysis.cov=analysis.cov, method.borrow=method.borrow,
#'   chains=1, iter=100, nsim=nsim)
#' @import overlapping stats
#' @export

psborrow.t2e <- function(
  n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, n.EC, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
  psmatch.cov,
  method.psest="glm", method.pslink="logit",
  method.whomatch, method.matching, method.psorder, n.boot=100,
  analysis.cov, method.borrow,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim)
{
  reject <- NULL
  theta  <- NULL
  ov     <- NULL

  for(ss in 1:nsim){

    indata <- trial.simulation.t2e(
      n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
      n.ECp=n.ECp, nevent.ECp=nevent.ECp, accrual=accrual,
      out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
      cov.C=cov.C, cov.cor.C=cov.cor.C,
      cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect)

    f1 <- stats::as.formula(paste("study~",paste(psmatch.cov,collapse="+"),sep=""))

    out.psmatch <- psmatch(
      formula=f1, data=indata, n.EC=n.EC,
      method.psest=method.psest, method.pslink=method.pslink,
      method.whomatch=method.whomatch, method.matching=method.matching,
      method.psorder=method.psorder, n.boot=n.boot)

    indata.match <- rbind(indata[indata$study==1,],indata[out.psmatch$subjid.EC,])

    f2 <- stats::as.formula(paste("Surv(time,status)~",paste(analysis.cov,collapse="+"),sep=""))

    out.commensurate <- commensurate.t2e(
      formula=f2, data=indata.match, method.borrow=method.borrow,
      chains=chains, iter=iter, warmup=warmup, thin=thin,
      alternative=alternative, sig.level=sig.level)

    reject <- rbind(reject,data.frame(sim=ss,out.commensurate$reject))
    theta  <- rbind(theta,data.frame(sim=ss,out.commensurate$theta))

    data.ps <- out.psmatch$data.ps

    Z1 <- data.ps[(data.ps$study==1)&(data.ps$treat==1),]
    Z2 <- rbind(data.ps[(data.ps$study==1)&(data.ps$treat==0),],data.ps[out.psmatch$subjid.EC,])

    pslist <- list(Z1=Z1[,"ps"],Z2=Z2[,"ps"])
    ov.ps  <- (overlap(pslist,boundaries=c(0,1)))$OV
    ov     <- rbind(ov,data.frame(sim=ss,factor="ps",type="continuous",comparison="CTvsCC.EC",ov=ov.ps))

    for(i in 1:length(psmatch.cov)){
      covlist  <- list(Z1=Z1[,psmatch.cov[i]],Z2=Z2[,psmatch.cov[i]])
      if(cov.C[[i]]$dist=="norm"){
        ov.cov <- (overlap(covlist))$OV
        ov     <- rbind(ov,data.frame(sim=ss,factor=psmatch.cov[i],type="continuous",comparison="CTvsCC.EC",ov=ov.cov))
      }else if(cov.C[[i]]$dist=="binom"){
        ov.cov <- mean(covlist$Z1)-mean(covlist$Z2)
        ov     <- rbind(ov,data.frame(sim=ss,factor=psmatch.cov[i],type="binary",comparison="CTvsCC.EC",ov=ov.cov))
      }
    }

    Z1 <- data.ps[(data.ps$study==1)&(data.ps$treat==0),]
    Z2 <- data.ps[out.psmatch$subjid.EC,]

    pslist <- list(Z1=Z1[,"ps"],Z2=Z2[,"ps"])
    ov.ps  <- (overlap(pslist,boundaries=c(0,1)))$OV
    ov     <- rbind(ov,data.frame(sim=ss,factor="ps",type="continuous",comparison="CCvsEC",ov=ov.ps))

    for(i in 1:length(psmatch.cov)){
      covlist  <- list(Z1=Z1[,psmatch.cov[i]],Z2=Z2[,psmatch.cov[i]])
      if(cov.C[[i]]$dist=="norm"){
        ov.cov <- (overlap(covlist))$OV
        ov     <- rbind(ov,data.frame(sim=ss,factor=psmatch.cov[i],type="continuous",comparison="CCvsEC",ov=ov.cov))
      }else if(cov.C[[i]]$dist=="binom"){
        ov.cov <- mean(covlist$Z1)-mean(covlist$Z2)
        ov     <- rbind(ov,data.frame(sim=ss,factor=psmatch.cov[i],type="binary",comparison="CCvsEC",ov=ov.cov))
      }
    }
  }

  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  t.theta <- log(out.lambda.CT/out.lambda.CC)

  return(list(reject=reject,theta=theta,ov=ov,
              n.CT=n.CT,n.CC=n.CC,n.ECp=n.ECp,n.EC=n.EC,drift=driftHR,
              true.theta=t.theta,
              method.psest=method.psest,method.pslink=method.pslink,
              method.whomatch=method.whomatch,
              method.matching=method.matching,method.psorder=method.psorder))
}

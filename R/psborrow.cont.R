
#' Simulation study of hybrid control design with Bayesian dynamic borrowing
#' incorporating propensity score matched external control: continuous outcome
#'
#' Bayesian dynamic borrowing is implemented for clinical trials with hybrid
#' control design, where the concurrent control is augmented by external
#' control. The external control is selected from external control pool using
#' propensity score matching. Commensurate power prior is used for Bayesian
#' dynamic borrowing, and half-normal or half-Cauchy commensurate prior is
#' available. No borrowing and full borrowing can also be implemented.
#' The continuous outcome is applicable.
#' @usage
#' psborrow.cont(
#'   n.CT, n.CC, n.ECp, n.EC,
#'   out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
#'   cov.C, cov.cor.C, cov.effect.C,
#'   cov.EC, cov.cor.EC, cov.effect.EC,
#'   psmatch.cov,
#'   method.psest="glm", method.pslink="logit",
#'   method.whomatch, method.matching, method.psorder,
#'   analysis.cov, method.borrow,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025,
#'   nsim, seed=sample.int(.Machine$integer.max,1))
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param n.ECp Number of patients in external control pool.
#' @param n.EC Number of patients in external control.
#' @param out.mean.CT True mean of outcome in concurrent treatment.
#' @param out.sd.CT True sd of outcome in concurrent treatment.
#' @param out.mean.CC True mean of outcome in concurrent control.
#' @param out.sd.CC True as of outcome in concurrent control.
#' @param driftdiff Difference between external control and concurrent control
#' for which the bias should be plotted.
#' @param out.sd.EC True as of outcome in external control.
#' @param cov.C List of covariate distributions for concurrent treatment and
#' control.
#' @param cov.cor.C Matrix of correlation coefficients for each pair of
#' covariate for concurrent treatment and control, specified as Gaussian copula
#' parameter.
#' @param cov.effect.C Vector of covariate effects for concurrent treatment and
#' control, specified as hazard ratio.
#' @param cov.EC List of covariate distributions for external control.
#' @param cov.cor.EC Matrix of correlation coefficients for each pair of
#' covariate for external control, specified as Gaussian copula parameter.
#' @param cov.effect.EC Vector of covariate effects for external control
#' control, specified as hazard ratio.
#' @param psmatch.cov psmatch.cov.
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
#' \code{conc.contl} matching concurrent control to external control pool;
#' \code{conc.treat} matching concurrent treatment to external control pool;
#' \code{conc.all} matching concurrent treatment plus concurrent control to
#' external control pool; \code{treat2contl} matching concurrent treatment to
#' concurrent control plus external control pool.
#' @param method.matching Matching method. Allowable options include
#' \code{"optimal"} for optimal matching; \code{"nearest"} for nearest neighbor
#' matching without replacement.
#' @param method.psorder Order that the matching takes place. Allowable options
#' include \code{"largest"}, where matching takes place in descending order of
#' propensity score; \code{"smallest"}, where matching takes place in ascending
#' order of propensity score; \code{"random"}, where matching takes place in a
#' random order; \code{"data"}, where matching takes place based on the order
#' of units in the data.
#' @param n.boot n.boot.
#' @param analysis.cov analysis.cov.
#' @param method.borrow List of information borrowing method. \code{"noborrow"}
#' uses the concurrent data only. \code{"fullborrow"} uses the external control
#' data without discounting. \code{"cauchy"} uses the commensurate prior to
#' dynamically borrow the external control data, and the commensurate parameter
#' is assumed to follow half-Cauchy distribution. \code{"normal"} uses the
#' commensurate prior to dynamically borrow the external control data, and the
#' commensurate parameter is assumed to follow half-normal distribution.
#' \code{"cauchy"} and \code{"normal"} require to specify the scale parameter of
#' half-Cauchy and half-normal distribution respectively.
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
#' @param seed Random seed.
#' @return
#' \item{reject}{\code{TRUE} when significant; otherwise \code{FALSE}.}
#' \item{theta}{Posterior mean, median, and sd of log hazard ratio.}
#' \item{ov.ps}{Overlapping coefficient of propensity score between concurrent
#' treatment and concurrent/external control}
#' \item{ov.cov}{Overlapping coefficient of continuous covariate and difference
#' of proportion of binary covariate between concurrent treatment and
#' concurrent/external control}
#' @examples
#' n.CT       <- 100
#' n.CC       <- 50
#' n.ECp      <- 1000
#' n.EC       <- 50
#'
#' out.mean.CT <- 0
#' out.sd.CT   <- 1
#' out.mean.CC <- 0
#' out.sd.CC   <- 1
#' driftdiff   <- 0
#' out.sd.EC   <- 1
#'
#' cov.C <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'               list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.C <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.C <- c(0.1,0.1)
#'
#' cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'                list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.EC <- rbind(c(  1,0.1),
#'                     c(0.1,  1))
#'
#' cov.effect.EC <- c(0.1,0.1)
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
#' psborrow.cont(
#'   n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp, n.EC=n.EC,
#'   out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
#'   out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
#'   driftdiff=driftdiff, out.sd.EC=out.sd.EC,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC,
#'   psmatch.cov=psmatch.cov, method.whomatch=method.whomatch,
#'   method.matching=method.matching, method.psorder=method.psorder,
#'   analysis.cov=analysis.cov, method.borrow=method.borrow,
#'   nsim=nsim, seed=100)
#' @import overlapping
#' @export

psborrow.cont <- function(
  n.CT, n.CC, n.ECp, n.EC,
  out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
  cov.C, cov.cor.C, cov.effect.C,
  cov.EC, cov.cor.EC, cov.effect.EC,
  psmatch.cov,
  method.psest="glm", method.pslink="logit",
  method.whomatch, method.matching, method.psorder, n.boot=100,
  analysis.cov, method.borrow,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025,
  nsim, seed=sample.int(.Machine$integer.max,1))
{
  reject <- NULL
  theta  <- NULL
  ov     <- NULL

  set.seed(seed)

  for(ss in 1:nsim){

    indata <- trial.simulation.cont(
      n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
      out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
      out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
      driftdiff=driftdiff, out.sd.EC=out.sd.EC,
      cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
      cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC)

    f1 <- as.formula(paste("study~",paste(psmatch.cov,collapse="+"),sep=""))

    out.psmatch <- psmatch(
      formula=f1, data=indata, n.EC=n.EC,
      method.psest=method.psest, method.pslink=method.pslink,
      method.whomatch=method.whomatch, method.matching=method.matching,
      method.psorder=method.psorder)

    indata.match <- rbind(indata[indata$study==1,],indata[out.psmatch$subjid.EC,])

    f2 <- as.formula(paste("y~",paste(analysis.cov,collapse="+"),sep=""))

    out.commensurate <- commensurate.cont(
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

  t.theta <- out.mean.CT-out.mean.CC

  return(list(reject=reject,theta=theta,ov=ov,
              n.CT=n.CT,n.CC=n.CC,n.ECp=n.ECp,n.EC=n.EC,drift=driftdiff,
              true.theta=t.theta,
              method.psest=method.psest,method.pslink=method.pslink,
              method.whomatch=method.whomatch,
              method.matching=method.matching,method.psorder=method.psorder))
}

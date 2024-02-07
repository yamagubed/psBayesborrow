
#' Summarizing simulation study results
#'
#' Simulation study results are summarized.
#' @usage
#' psborrow.summary(object)
#' @param object List of simulation results.
#' @return
#' The \code{psborrow.summary} returns a list containing the following objects:
#' \item{res.out}{Data frame containing (1) rate that the null hypothesis is
#' rejected in Bayesian one-sided hypothesis testing, (2) bias in posterior mean
#' of treatment effect (log hazard ratio, log odds ratio, or mean difference),
#' (3) empirical standard deviation (\code{EmpSD}), derived by taking the
#' standard deviation of posterior mean of treatment effect that each simulated
#' trial yielded, (4) model standard deviation ((\code{ModSD})), derived by
#' taking the average of posterior standard deviation of treatment effect that
#' each simulated trial yielded, (5) ratio of \code{ModSD} to \code{EmpSD},
#' (6) Coverage probability of credible interval.}
#' \item{res.cov}{Data frame containing (1) mean, median, and standard deviation
#' of overlapping coefficients of propensity score densities between treatment
#' versus concurrent control plus external control and between concurrent
#' control versus external control, (2) mean, median, and standard deviation
#' of overlapping coefficients of continuous covariate densities between treatment
#' versus concurrent control plus external control and between concurrent
#' control versus external control, and (3) mean, median, and standard deviation
#' of rate differences of binary covariate between treatment versus concurrent
#' control plus external control and between concurrent control versus external
#' control.}
#' \item{n.CT}{Number of patients in treatment group in the current trial.}
#' \item{n.CC}{Number of patients in concurrent control group in the current
#' trial.}
#' \item{n.ECp}{Number of patients in external control pool.}
#' \item{n.EC}{Number of patients in external control.}
#' \item{drift}{Hazard ratio, odds ratio, or mean difference between concurrent
#' and external control for which the bias should be plotted. The measure
#' depends on the outcome type.}
#' \item{method.psest}{Method of estimating the propensity score.}
#' \item{method.pslink}{Link function used in estimating the propensity score.}
#' \item{method.whomatch}{Options of who to match.}
#' \item{method.matching}{Matching method.}
#' \item{method.psorder}{Order that the matching takes place when a nearest
#' neighbor matching is used.}
#' @import stats
#' @export
psborrow.summary <- function(object)
{
  reject  <- object$reject
  theta   <- object$theta
  colflg  <- colnames(object$reject)[-c(1,2)]
  t.theta <- object$true.theta

  rrate    <- apply(reject[,colflg],2,mean)
  bias     <- apply(theta[theta$measure=="mean",colflg]-t.theta,2,mean)
  empsd    <- apply(theta[theta$measure=="mean",colflg],2,stats::sd)
  modsd    <- apply(theta[theta$measure=="sd",colflg],2,mean)
  modempsd <- modsd/empsd
  covprob  <- apply((theta[theta$measure=="lcri",colflg]<t.theta)*(theta[theta$measure=="ucri",colflg]>t.theta),2,mean)

  res.out <- data.frame(measure=c("Reject","Bias","EmpSD","ModSD","ModSD/EmpSD","Coverage Prob."),
                        round(rbind(rrate,bias,empsd,modsd,modempsd,covprob),digits=3))

  colnames(res.out) <- c("measure",colflg)
  rownames(res.out) <- NULL

  ov <- object$ov

  ov.ft   <- unique(ov[,c("factor","type","comparison")])
  l.ov.ft <- nrow(ov.ft)
  res.cov <- NULL
  for(i in 1:l.ov.ft){
    anaflg <- (ov$factor==ov.ft$factor[i])&(ov$type==ov.ft$type[i])&(ov$comparison==ov.ft$comparison[i])

    res.cov <- rbind(res.cov,data.frame(
                       factor     = ov.ft$factor[i],
                       type       = ov.ft$type[i],
                       comparison = ov.ft$comparison[i],
                       mean       = round(mean(ov[anaflg,"ov"]),digits=3),
                       median     = round(stats::median(ov[anaflg,"ov"]),digits=3),
                       sd         = round(stats::sd(ov[anaflg,"ov"]),digits=3)
                    ))
  }

  result <- list(res.out         = res.out,
                 res.cov         = res.cov,
                 n.CT            = object$n.CT,
                 n.CC            = object$n.CC,
                 n.ECp           = object$n.ECp,
                 n.EC            = object$n.EC,
                 drift           = object$drift,
                 method.psest    = object$method.psest,
                 method.pslink   = object$method.pslink,
                 method.whomatch = object$method.whomatch,
                 method.matching = object$method.matching,
                 method.psorder  = object$method.psorder)

  class(result) <- "psborrow.summary"
  result
}

#' @export
print.psborrow.summary <- function(x,...)
{
  ptab1  <- as.matrix(x$res.out[,-1])
  rnames <- c("Reject","Bias","EmpSD","ModSD","EmpSD/ModSD","Coverage Prob.")
  dimnames(ptab1) <- list(rnames,names(x$res.out[,-1]))

  ptab2  <- as.matrix(x$res.cov[,-1])
  rnames <- x$res.cov[,1]
  dimnames(ptab2) <- list(rnames,names(x$res.cov[,-1]))

  ptab3 <- rbind(
    x$n.CT,
    x$n.CC,
    x$n.ECp,
    x$n.EC,
    x$drift)
  rnames <- c("No. of conc. treatment",
              "No. of conc. control",
              "No. of ext. control pool",
              "No. of ext. control",
              "Prior data conflict")
  dimnames(ptab3) <- list(rnames,"")

  if(x$method.matching=="optimal"){
    ptab4 <- rbind(
      x$method.psest,
      x$method.pslink,
      x$method.whomatch,
      x$method.matching)
    rnames <- c("PS estimation",
                "Link function",
                "Who to match",
                "Matching method")
    dimnames(ptab4) <- list(rnames,"")
  }else{
    ptab4 <- rbind(
      x$method.psest,
      x$method.pslink,
      x$method.whomatch,
      x$method.matching,
      x$method.psorder)
    rnames <- c("PS estimation",
                "Link function",
                "Who to match",
                "Matching method",
                "Matching order")
    dimnames(ptab4) <- list(rnames,"")
  }

  cat("\n")
  cat("Treatment effect estimation:\n")
  print.default(ptab1,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Overlapping coefficient:\n")
  print.default(ptab2,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Trial design settings:\n")
  print.default(ptab3,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Method:\n")
  print.default(ptab4,print.gap=2L,quote=FALSE)
  cat("\n")
}

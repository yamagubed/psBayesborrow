
#' Simulating time-to-event data for current trial and external control
#'
#' A two-arm randomized clinical trial with a time-to-event outcome, which is
#' augmented by external control data, is simulated.
#' @usage
#' trial.simulation.t2e(
#'  n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, accrual,
#'  out.mevent.CT, out.mevent.CC, driftHR,
#'  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
#'  seed=sample.int(.Machine$integer.max,1))
#' @param n.CT Number of patients in treatment group in the current trial.
#' @param n.CC Number of patients in concurrent control group in the current
#' trial.
#' @param nevent.C Number of events in treatment and concurrent control group
#' in the current trial.
#' @param n.ECp Number of patients in external control pool.
#' @param nevent.ECp Number of events in external control pool.
#' @param accrual Accrual rate, defined as the number of enrolled patients per
#' month.
#' @param out.mevent.CT True median time to event in treatment group in the
#' current trial.
#' @param out.mevent.CC True median time to event in concurrent control group
#' in the current trial.
#' @param driftHR Hazard ratio between concurrent and external control for
#' which the bias should be plotted (hazard in external control divided by
#' hazard in concurrent control).
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
#' @param cov.effect Vector of covariate effects on the outcome , specified as
#' hazard ratio per one unit increase in continuous covariates or as hazard
#' ratio between categories for binary covariates.
#' @param seed Setting a seed.
#' @details The time to event outcome is assumed to follow a Weibull
#' distribution. Given more than one covariates with their effects on the
#' outcome, a Weibull proportional hazards model is constructed for data
#' generation. The data frame generated include the time-to-event outcome data
#' and covariates for \code{n.CT} and \code{n.CC} patients in treatment and
#' concurrent control group in the current trial respectively, and \code{n.ECp}
#' patients in external control pool. One record per patient. More than one
#' covariates must be specified.
#' @return
#' The \code{trial.simulation.t2e} returns a data frame containing the
#' following variables:
#' \item{study}{Study indicator (0 for external control, and 1 for current
#' trial)}
#' \item{treat}{Treatment indicator (0 for concurrent and external control, and
#' 1 for treatment)}
#' \item{time}{Time to event or censoring}
#' \item{status}{Censoring (0 for censored, and 1 for event occurred)}
#' \item{column name specified}{Covariate of interest}
#' @examples
#' n.CT       <- 100
#' n.CC       <- 50
#' nevent.C   <- 100
#' n.ECp      <- 1000
#' nevent.ECp <- 800
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
#' cov.effect <- c(0.8,0.8)
#'
#' trial.simulation.t2e(
#'    n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'    n.ECp=n.ECp, nevent.ECp=nevent.ECp, accrual=accrual,
#'    out.mevent.CT, out.mevent.CC, driftHR,
#'    cov.C=cov.C, cov.cor.C=cov.cor.C,
#'    cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect, seed=100)
#' @import stats
#' @export

trial.simulation.t2e <- function(
  n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
  seed=sample.int(.Machine$integer.max,1))
{
  set.seed(seed)

  ncov          <- length(cov.C)
  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC

  lcov.effect  <- (-log(cov.effect))

  marg.C  <- NULL
  marg.EC <- NULL
  mean.C  <- NULL
  mean.EC <- NULL
  cov.lab <- NULL

  for(i in 1:ncov){
    if(cov.C[[i]]$dist=="norm"){
      marg.C  <- append(marg.C, list(list(dist=cov.C[[i]]$dist, parm=list(mean=cov.C[[i]]$mean, sd=cov.C[[i]]$sd))))
      marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(mean=cov.EC[[i]]$mean,sd=cov.EC[[i]]$sd))))

      mean.C  <- c(mean.C, cov.C[[i]]$mean)
      mean.EC <- c(mean.EC,cov.EC[[i]]$mean)

      cov.lab <- c(cov.lab,cov.C[[i]]$lab)
    }else if(cov.C[[i]]$dist=="binom"){
      marg.C  <- append(marg.C, list(list(dist=cov.C[[i]]$dist, parm=list(size=1,prob=cov.C[[i]]$prob))))
      marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(size=1,prob=cov.EC[[i]]$prob))))

      mean.C  <- c(mean.C, cov.C[[i]]$prob)
      mean.EC <- c(mean.EC,cov.EC[[i]]$prob)

      cov.lab <- c(cov.lab,cov.C[[i]]$lab)
    }
  }

  int.C   <- log(1/out.lambda.CC)-sum(mean.C*lcov.effect)
  t.theta <- log(out.lambda.CC/out.lambda.CT)

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  sigma.CT  <- exp(int.C+t.theta             +apply(data.cov.CT, 1,function(x){sum(x*lcov.effect)}))
  sigma.CC  <- exp(int.C                     +apply(data.cov.CC, 1,function(x){sum(x*lcov.effect)}))
  sigma.ECp <- exp(int.C        -log(driftHR)+apply(data.cov.ECp,1,function(x){sum(x*lcov.effect)}))

  data.CT  <- cbind(stats::rweibull(n.CT, shape=1,scale=sigma.CT), data.cov.CT)
  data.CC  <- cbind(stats::rweibull(n.CC, shape=1,scale=sigma.CC), data.cov.CC)
  data.ECp <- cbind(stats::rweibull(n.ECp,shape=1,scale=sigma.ECp),data.cov.ECp)

  enroll.period.C <- floor((n.CT+n.CC)/accrual)
  mn.CT <- floor(n.CT/enroll.period.C)
  mn.CC <- floor(n.CC/enroll.period.C)

  enroll.time.CT <- NULL
  enroll.time.CC <- NULL
  for(i in 1:enroll.period.C){
    e.st <- i-1
    e.en <- i
    enroll.time.CT <- c(enroll.time.CT,stats::runif(mn.CT,e.st,e.en))
    enroll.time.CC <- c(enroll.time.CC,stats::runif(mn.CC,e.st,e.en))
  }
  enroll.time.CT <- c(enroll.time.CT,stats::runif(n.CT-length(enroll.time.CT),enroll.period.C,enroll.period.C+1))
  enroll.time.CC <- c(enroll.time.CC,stats::runif(n.CC-length(enroll.time.CC),enroll.period.C,enroll.period.C+1))

  enroll.period.ECp <- floor(n.ECp/accrual)

  enroll.time.ECp <- NULL
  for(i in 1:enroll.period.ECp){
    e.st <- i-1
    e.en <- i
    enroll.time.ECp <- c(enroll.time.ECp,stats::runif(accrual,e.st,e.en))
  }
  enroll.time.ECp <- c(enroll.time.ECp,stats::runif(n.ECp-length(enroll.time.ECp),enroll.period.ECp,enroll.period.ECp+1))

  obs.time.CT  <- data.CT[,1] +enroll.time.CT
  obs.time.CC  <- data.CC[,1] +enroll.time.CC
  obs.time.ECp <- data.ECp[,1]+enroll.time.ECp

  last.sub.C   <- sort(c(obs.time.CT,obs.time.CC))[nevent.C]
  last.sub.ECp <- sort(obs.time.ECp)[nevent.ECp]

  censor.CT  <- as.numeric(obs.time.CT>last.sub.C)
  censor.CC  <- as.numeric(obs.time.CC>last.sub.C)
  censor.ECp <- as.numeric(obs.time.ECp>last.sub.ECp)

  cenf <- (obs.time.CT>last.sub.C)
  data.CT[cenf,1] <- data.CT[cenf,1]-(obs.time.CT[cenf]-last.sub.C)

  cenf <- (obs.time.CC>last.sub.C)
  data.CC[cenf,1] <- data.CC[cenf,1]-(obs.time.CC[cenf]-last.sub.C)

  cenf <- (obs.time.ECp>last.sub.ECp)
  data.ECp[cenf,1] <- data.ECp[cenf,1]-(obs.time.ECp[cenf]-last.sub.ECp)

  censor.CT <- censor.CT[data.CT[,1]>0]
  data.CT   <- data.CT[data.CT[,1]>0,]

  censor.CC <- censor.CC[data.CC[,1]>0]
  data.CC   <- data.CC[data.CC[,1]>0,]

  censor.ECp <- censor.ECp[data.ECp[,1]>0]
  data.ECp   <- data.ECp[data.ECp[,1]>0,]

  outdata <- rbind(
    data.frame(study=1,treat=1,time=data.CT[,1], status=1-censor.CT, data.CT[,-1,drop=FALSE]),
    data.frame(study=1,treat=0,time=data.CC[,1], status=1-censor.CC, data.CC[,-1,drop=FALSE]),
    data.frame(study=0,treat=0,time=data.ECp[,1],status=1-censor.ECp,data.ECp[,-1,drop=FALSE]))

  colnames(outdata) <- c("study","treat","time","status",cov.lab)

  return(outdata)
}


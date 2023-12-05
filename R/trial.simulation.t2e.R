
#' Simulating a trial with time-to-event outcome
#'
#' A trial with time-to-event outcome is simulated. The data returned include
#' the time to event of interest, an indicator of censoring, and covariate for
#' patients in concurrent treatment and control, and external control pool.
#' @usage
#' trial.simulation.t2e(
#'  n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, accrual,
#'  out.mevent.CT, out.mevent.CC, driftHR,
#'  cov.C, cov.cor.C, cov.effect.C,
#'  cov.EC, cov.cor.EC, cov.effect.EC)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param nevent.C Number of events in concurrent treatment and control.
#' @param n.ECp Number of patients in external control pool.
#' @param nevent.ECp Number of events in external control pool.
#' @param accrual Accrual rate (number of enrolled patients per month).
#' @param out.mevent.CT True median time to event in concurrent treatment.
#' @param out.mevent.CC True median time to event in concurrent control.
#' @param driftHR Hazard ratio between concurrent and external control for
#' which the bias should be plotted.
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
#' @return
#' \item{outdata}{Dataset of a simulated trial.}
#' \item{ncov}{Number of covariate.}
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
#' trial.simulation.t2e(
#'    n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'    n.ECp=n.ECp, nevent.ECp=nevent.ECp, accrual=accrual,
#'    out.mevent.CT, out.mevent.CC, driftHR,
#'    cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'    cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC)
#' @export

trial.simulation.t2e <- function(
  n.CT, n.CC, nevent.C, n.ECp, nevent.ECp, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C, cov.cor.C, cov.effect.C,
  cov.EC, cov.cor.EC, cov.effect.EC)
{
  ncov          <- length(cov.C)
  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  out.lambda.EC <- out.lambda.CC*driftHR

  lcov.effect.C  <- (-log(cov.effect.C))
  lcov.effect.EC <- (-log(cov.effect.EC))

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

  int.C   <- log(1/out.lambda.CC)-sum(mean.C*lcov.effect.C)
  int.EC  <- log(1/out.lambda.EC)-sum(mean.EC*lcov.effect.EC)
  t.theta <- log(out.lambda.CC/out.lambda.CT)

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  sigma.CT  <- exp(int.C +t.theta+apply(data.cov.CT, 1,function(x){sum(x*lcov.effect.C)}))
  sigma.CC  <- exp(int.C         +apply(data.cov.CC, 1,function(x){sum(x*lcov.effect.C)}))
  sigma.ECp <- exp(int.EC        +apply(data.cov.ECp,1,function(x){sum(x*lcov.effect.EC)}))

  data.CT  <- cbind(rweibull(n.CT, shape=1,scale=sigma.CT), data.cov.CT)
  data.CC  <- cbind(rweibull(n.CC, shape=1,scale=sigma.CC), data.cov.CC)
  data.ECp <- cbind(rweibull(n.ECp,shape=1,scale=sigma.ECp),data.cov.ECp)

  enroll.period.C <- floor((n.CT+n.CC)/accrual)
  mn.CT <- floor(n.CT/enroll.period.C)
  mn.CC <- floor(n.CC/enroll.period.C)

  enroll.time.CT <- NULL
  enroll.time.CC <- NULL
  for(i in 1:enroll.period.C){
    e.st <- i-1
    e.en <- i
    enroll.time.CT <- c(enroll.time.CT,runif(mn.CT,e.st,e.en))
    enroll.time.CC <- c(enroll.time.CC,runif(mn.CC,e.st,e.en))
  }
  enroll.time.CT <- c(enroll.time.CT,runif(n.CT-length(enroll.time.CT),enroll.period.C,enroll.period.C+1))
  enroll.time.CC <- c(enroll.time.CC,runif(n.CC-length(enroll.time.CC),enroll.period.C,enroll.period.C+1))

  enroll.period.ECp <- floor(n.ECp/accrual)

  enroll.time.ECp <- NULL
  for(i in 1:enroll.period.ECp){
    e.st <- i-1
    e.en <- i
    enroll.time.ECp <- c(enroll.time.ECp,runif(accrual,e.st,e.en))
    }
  enroll.time.ECp <- c(enroll.time.ECp,runif(n.ECp-length(enroll.time.ECp),enroll.period.ECp,enroll.period.ECp+1))

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
    data.frame(study=1,treat=1,time=data.CT[,1], status=1-censor.CT, data.CT[,-1]),
    data.frame(study=1,treat=0,time=data.CC[,1], status=1-censor.CC, data.CC[,-1]),
    data.frame(study=0,treat=0,time=data.ECp[,1],status=1-censor.ECp,data.ECp[,-1]))

  colnames(outdata) <- c("study","treat","time","status",cov.lab)

  return(outdata)
}


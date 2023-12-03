
#' Simulating a trial with binary outcome
#'
#' A trial with binary outcome is simulated. The data returned include
#' the time to event of interest, an indicator of censoring, and covariate for
#' patients in concurrent treatment and control, and external control pool.
#' @usage
#' trial.simulation.bin(
#'  n.CT, n.CC, n.ECp,
#'  out.prob.CT, out.prob.CC, driftOR,
#'  cov.C, cov.cor.C, cov.effect.C,
#'  cov.EC, cov.cor.EC, cov.effect.EC)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param n.ECp Number of patients in external control pool.
#' @param out.prob.CT True probability of outcome in concurrent treatment.
#' @param out.prob.CC True probability of outcome in concurrent control.
#' @param driftOR OR between external control and concurrent control for which
#' the bias should be plotted.
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
#' n.ECp      <- 1000
#'
#' out.prob.CT <- 0.2
#' out.prob.CC <- 0.2
#' driftOR     <- 1.0
#'
#' cov.C <- list(list(dist="norm",mean=0,sd=1),
#'               list(dist="binom",prob=0.4))
#'
#' cov.cor.C <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.C <- c(0.1,0.1)
#'
#' cov.EC <- list(list(dist="norm",mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cov.cor.EC <- rbind(c(  1,0.1),
#'                     c(0.1,  1))
#'
#' cov.effect.EC <- c(0.1,0.1)
#'
#' trial.simulation.bin(
#'    n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
#'    out.prob.CT=out.prob.CT, out.prob.CC=out.prob.CC, driftOR=driftOR,
#'    cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'    cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC)
#' @import boot
#' @export

trial.simulation.bin <- function(
  n.CT, n.CC, n.ECp,
  out.prob.CT, out.prob.CC, driftOR,
  cov.C, cov.cor.C, cov.effect.C,
  cov.EC, cov.cor.EC, cov.effect.EC)
{
  ncov        <- length(cov.C)
  out.prob.EC <- (driftOR*out.prob.CC/(1-out.prob.CC))/(1+driftOR*out.prob.CC/(1-out.prob.CC))

  lcov.effect.C  <- log(cov.effect.C)
  lcov.effect.EC <- log(cov.effect.EC)

  marg.C  <- NULL
  marg.EC <- NULL
  mean.C  <- NULL
  mean.EC <- NULL

  for(i in 1:ncov){
    if(cov.C[[i]]$dist=="norm"){
      marg.C  <- append(marg.C, list(list(dist=cov.C[[i]]$dist, parm=list(mean=cov.C[[i]]$mean, sd=cov.C[[i]]$sd))))
      marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(mean=cov.EC[[i]]$mean,sd=cov.EC[[i]]$sd))))

      mean.C  <- c(mean.C, cov.C[[i]]$mean)
      mean.EC <- c(mean.EC,cov.EC[[i]]$mean)
    }else if(cov.C[[i]]$dist=="binom"){
      marg.C  <- append(marg.C, list(list(dist=cov.C[[i]]$dist, parm=list(size=1,prob=cov.C[[i]]$prob))))
      marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(size=1,prob=cov.EC[[i]]$prob))))

      mean.C  <- c(mean.C, cov.C[[i]]$prob)
      mean.EC <- c(mean.EC,cov.EC[[i]]$prob)
    }
  }

  int.C   <- log(out.prob.CC/(1-out.prob.CC))-sum(mean.C*lcov.effect.C)
  int.EC  <- log(out.prob.EC/(1-out.prob.EC))-sum(mean.EC*lcov.effect.EC)
  t.theta <- log(out.prob.CT/(1-out.prob.CT))-log(out.prob.CC/(1-out.prob.CC))

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  p.CT  <- boot::inv.logit(int.C +t.theta+apply(data.cov.CT, 1,function(x){sum(x*lcov.effect.C)}))
  p.CC  <- boot::inv.logit(int.C         +apply(data.cov.CC, 1,function(x){sum(x*lcov.effect.C)}))
  p.ECp <- boot::inv.logit(int.EC        +apply(data.cov.ECp,1,function(x){sum(x*lcov.effect.EC)}))

  data.CT  <- cbind(rbinom(n.CT, 1,p.CT), data.cov.CT)
  data.CC  <- cbind(rbinom(n.CC, 1,p.CC), data.cov.CC)
  data.ECp <- cbind(rbinom(n.ECp,1,p.ECp),data.cov.ECp)

  outdata <- rbind(
    data.frame(study=1,treat=1,y=data.CT[,1], data.CT[,-1]),
    data.frame(study=1,treat=0,y=data.CC[,1], data.CC[,-1]),
    data.frame(study=0,treat=0,y=data.ECp[,1],data.ECp[,-1]))

  return(outdata)
}


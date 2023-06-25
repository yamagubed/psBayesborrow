
#' Simulating a trial with continuous outcome
#'
#' A trial with continuous outcome is simulated. The data returned include
#' the time to event of interest, an indicator of censoring, and covariate for
#' patients in concurrent treatment and control, and external control pool.
#' @usage
#' trial.simulation.cont(
#'  n.CT, n.CC, n.ECp,
#'  out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
#'  cov.C, cov.cor.C, cov.effect.C,
#'  cov.EC, cov.cor.EC, cov.effect.EC)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param n.ECp Number of patients in external control pool.
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
#' @return
#' \item{outdata}{Dataset of a simulated trial.}
#' \item{ncov}{Number of covariate.}
#' @examples
#' n.CT       <- 100
#' n.CC       <- 50
#' n.ECp      <- 1000
#'
#' out.mean.CT <- 0
#' out.sd.CT   <- 1
#' out.mean.CC <- 0
#' out.sd.CC   <- 1
#' driftdiff   <- 0
#' out.sd.EC   <- 1
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
#' trial.simulation.cont(
#'   n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
#'   out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
#'   out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
#'   driftdiff=driftdiff, out.sd.EC=out.sd.EC,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC)
#' @export

trial.simulation.cont <- function(
    n.CT, n.CC, n.ECp,
    out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
    cov.C, cov.cor.C, cov.effect.C,
    cov.EC, cov.cor.EC, cov.effect.EC)
{
  ncov        <- length(cov.C)
  out.mean.EC <- driftdiff+out.mean.CC

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

  int.C   <- out.mean.CC-sum(mean.C*cov.effect.C)
  int.EC  <- out.mean.EC-sum(mean.EC*cov.effect.EC)
  t.theta <- out.mean.CT-out.mean.CC

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  mu.CT  <- int.C +t.theta+apply(data.cov.CT, 1,function(x){sum(x*cov.effect.C)})
  mu.CC  <- int.C         +apply(data.cov.CC, 1,function(x){sum(x*cov.effect.C)})
  mu.ECp <- int.EC        +apply(data.cov.ECp,1,function(x){sum(x*cov.effect.EC)})

  data.CT  <- cbind(rnorm(n.CT, mean=mu.CT, sd=out.sd.CT),data.cov.CT)
  data.CC  <- cbind(rnorm(n.CC, mean=mu.CC, sd=out.sd.CC),data.cov.CC)
  data.ECp <- cbind(rnorm(n.ECp,mean=mu.ECp,sd=out.sd.EC),data.cov.ECp)

  outdata <- rbind(
    data.frame(study=1,treat=1,y=data.CT[,1], data.CT[,-1]),
    data.frame(study=1,treat=0,y=data.CC[,1], data.CC[,-1]),
    data.frame(study=0,treat=0,y=data.ECp[,1],data.ECp[,-1]))

  return(list(data=outdata,ncov=ncov))
}

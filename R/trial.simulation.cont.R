
#' Simulating continuous data for current trial and external control
#'
#' A two-arm randomized clinical trial with a continuous outcome, which is
#' augmented by external control data, is simulated.
#' @usage
#' trial.simulation.cont(
#'  n.CT, n.CC, n.ECp,
#'  out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
#'  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
#'  seed=sample.int(.Machine$integer.max,1))
#' @param n.CT Number of patients in treatment group in the current trial.
#' @param n.CC Number of patients in concurrent control group in the current
#' trial.
#' @param n.ECp Number of patients in external control pool.
#' @param out.mean.CT True mean of outcome in treatment group in the current
#' trial.
#' @param out.sd.CT True sd of outcome in treatment group in the current trial.
#' @param out.mean.CC True mean of outcome in concurrent control group in the
#' current trial.
#' @param out.sd.CC True sd of outcome in concurrent control group in the
#' current trial.
#' @param driftdiff Mean difference between concurrent and external control
#' for which the bias should be plotted (mean in external control minus mean
#' in concurrent control).
#' @param out.sd.EC True sd of outcome in external control.
#' @param cov.C List of covariate distributions for treatment and concurrent
#' control group in the current trial. Continuous and binary covariate are
#' applicable. The continuous covariate is assumed to follow a normal
#' distribution; for example specified as
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
#' mean change per one unit increase in continuous covariates or as mean change
#' between categories for binary covariates.
#' @param seed Setting a seed.
#' @details The continuous outcome is assumed to follow a normal distribution.
#' Given more than one covariates with their effects on the outcome, a normal
#' linear regression model is constructed for data generation. The data frame
#' generated include the continuous outcome data and covariates for \code{n.CT}
#' and \code{n.CC} patients in treatment and concurrent control group in the
#' current trial respectively, and \code{n.ECp} patients in external control
#' pool. One record per patient. More than one covariates must be specified.
#' @return
#' The \code{trial.simulation.cont} returns a data frame containing the
#' following variables:
#' \item{study}{Study indicator (0 for external control, and 1 for current
#' trial)}
#' \item{treat}{Treatment indicator (0 for concurrent and external control, and
#' 1 for treatment)}
#' \item{y}{Continuous outcome}
#' \item{column name specified}{Covariate of interest}
#' @examples
#' n.CT  <- 100
#' n.CC  <- 50
#' n.ECp <- 1000
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
#' cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'                list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.EC <- rbind(c(  1,0.1),
#'                     c(0.1,  1))
#'
#' cov.effect <- c(0.1,0.1)
#'
#' trial.simulation.cont(
#'   n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
#'   out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
#'   out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
#'   driftdiff=driftdiff, out.sd.EC=out.sd.EC,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C,
#'   cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect, seed=100)
#' @import stats
#' @export

trial.simulation.cont <- function(
    n.CT, n.CC, n.ECp,
    out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
    cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect,
    seed=sample.int(.Machine$integer.max,1))
{
  set.seed(seed)

  ncov <- length(cov.C)

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

  int.C   <- out.mean.CC-sum(mean.C*cov.effect)
  t.theta <- out.mean.CT-out.mean.CC

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  mu.CT  <- int.C+t.theta          +apply(data.cov.CT, 1,function(x){sum(x*cov.effect)})
  mu.CC  <- int.C                  +apply(data.cov.CC, 1,function(x){sum(x*cov.effect)})
  mu.ECp <- int.C        +driftdiff+apply(data.cov.ECp,1,function(x){sum(x*cov.effect)})

  data.CT  <- cbind(stats::rnorm(n.CT, mean=mu.CT, sd=out.sd.CT),data.cov.CT)
  data.CC  <- cbind(stats::rnorm(n.CC, mean=mu.CC, sd=out.sd.CC),data.cov.CC)
  data.ECp <- cbind(stats::rnorm(n.ECp,mean=mu.ECp,sd=out.sd.EC),data.cov.ECp)

  outdata <- rbind(
    data.frame(study=1,treat=1,y=data.CT[,1], data.CT[,-1,drop=FALSE]),
    data.frame(study=1,treat=0,y=data.CC[,1], data.CC[,-1,drop=FALSE]),
    data.frame(study=0,treat=0,y=data.ECp[,1],data.ECp[,-1,drop=FALSE]))

  colnames(outdata) <- c("study","treat","y",cov.lab)

  return(outdata)
}

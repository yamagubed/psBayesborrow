
#' Simulating binary data for current trial and external control
#'
#' A two-arm randomized clinical trial with a binary outcome, which is
#' augmented by external control data, is simulated.
#' @usage
#' trial.simulation.bin(
#'  n.CT, n.CC, n.ECp,
#'  out.prob.CT, out.prob.CC, driftOR,
#'  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect)
#' @param n.CT Number of patients in treatment group in the current trial.
#' @param n.CC Number of patients in concurrent control group in the current
#' trial.
#' @param n.ECp Number of patients in external control pool.
#' @param out.prob.CT True rate of outcome in treatment group in the current
#' trial.
#' @param out.prob.CC True rate of outcome in concurrent control group in the
#' current trial.
#' @param driftOR Odds ratio between concurrent and external control
#' for which the bias should be plotted.
#' @param cov.C List of covariate distributions for treatment and concurrent
#' control group in the current trial. Continuous and binary covariate are
#' applicable. The continuous covariate is assumed to follow a normal
#' distribution; for example, specified as
#' \code{list(dist="norm", mean=0, sd=1, lab="cov1")}. The binary covariate is
#' assumed to follow a binomial distribution; for example, specified as
#' \code{list(dist="binom", prob=0.4, lab="cov2")}. \code{lab} is the column
#' name of the covariate in the data frame generated.
#' @param cov.cor.C Matrix of correlation coefficients for each pair of
#' covariate for treatment and concurrent control in the current trial,
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
#' odds ratio per one unit increase in continuous covariates or as odds ratio
#' between categories for binary covariates.
#' @details The binary outcome is assumed to follow a binomial distribution.
#' Given more than one covariates with their effects on the outcome, a logistic
#' regression model is constructed for data generation. The data frame
#' generated include the binary outcome data and covariates for \code{n.CT}
#' and \code{n.CC} patients in treatment and concurrent control group in the
#' current trial respectively, and \code{n.ECp} patients in external control
#' pool. One record per patient. More than one covariates must be specified.
#' @return
#' The \code{trial.simulation.bin} returns a data frame containing the
#' following variables:
#' \item{study}{Study indicator (0 for external control, and 1 for current
#' trial)}
#' \item{treat}{Treatment indicator (0 for concurrent and external control, and
#' 1 for treatment)}
#' \item{y}{Binary outcome}
#' \item{column name specified}{Covariate of interest}
#' @examples
#' n.CT  <- 100
#' n.CC  <- 50
#' n.ECp <- 1000
#'
#' out.prob.CT <- 0.2
#' out.prob.CC <- 0.2
#' driftOR     <- 1.0
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
#' trial.simulation.bin(
#'    n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
#'    out.prob.CT=out.prob.CT, out.prob.CC=out.prob.CC, driftOR=driftOR,
#'    cov.C=cov.C, cov.cor.C=cov.cor.C,
#'    cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect)
#' @import boot stats
#' @export

trial.simulation.bin <- function(
  n.CT, n.CC, n.ECp,
  out.prob.CT, out.prob.CC, driftOR,
  cov.C, cov.cor.C, cov.EC, cov.cor.EC, cov.effect)
{
  ncov        <- length(cov.C)
  lcov.effect <- log(cov.effect)

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

  int.C   <- log(out.prob.CC/(1-out.prob.CC))-sum(mean.C*lcov.effect)
  t.theta <- log(out.prob.CT/(1-out.prob.CT))-log(out.prob.CC/(1-out.prob.CC))

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.EC <- cov.cor.EC[lower.tri(cov.cor.EC)]

  data.cov.CT  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CT)
  data.cov.CC  <- datagen(margdist=marg.C, corvec=cvec.C, nsim=n.CC)
  data.cov.ECp <- datagen(margdist=marg.EC,corvec=cvec.EC,nsim=n.ECp)

  p.CT  <- boot::inv.logit(int.C+t.theta             +apply(data.cov.CT, 1,function(x){sum(x*lcov.effect)}))
  p.CC  <- boot::inv.logit(int.C                     +apply(data.cov.CC, 1,function(x){sum(x*lcov.effect)}))
  p.ECp <- boot::inv.logit(int.C        +log(driftOR)+apply(data.cov.ECp,1,function(x){sum(x*lcov.effect)}))

  data.CT  <- cbind(stats::rbinom(n.CT, 1,p.CT), data.cov.CT)
  data.CC  <- cbind(stats::rbinom(n.CC, 1,p.CC), data.cov.CC)
  data.ECp <- cbind(stats::rbinom(n.ECp,1,p.ECp),data.cov.ECp)

  outdata <- rbind(
    data.frame(study=1,treat=1,y=data.CT[,1], data.CT[,-1]),
    data.frame(study=1,treat=0,y=data.CC[,1], data.CC[,-1]),
    data.frame(study=0,treat=0,y=data.ECp[,1],data.ECp[,-1]))

  colnames(outdata) <- c("study","treat","y",cov.lab)

  return(outdata)
}


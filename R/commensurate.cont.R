
#' Simulation study of hybrid control design with commensurate power prior:
#' continuous outcome
#'
#' Bayesian dynamic information borrowing with commensurate power prior is
#' implemented for hybrid control design, where the concurrent control is
#' augmented by historical control. No borrowing and full borrowing are also
#' implemented. The continuous outcome is applicable.
#' @usage
#' commensurate.cont(
#'   n.CT, n.CC, n.EC,
#'   out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
#'   cov.CT, cov.CC, cov.EC, cormat, method,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025, nsim=10)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param n.EC Number of patients in external control.
#' @param out.mean.CT True mean of outcome in concurrent treatment.
#' @param out.sd.CT True sd of outcome in concurrent treatment.
#' @param out.mean.CC True mean of outcome in concurrent control.
#' @param out.sd.CC True as of outcome in concurrent control.
#' @param driftdiff Difference between external control and concurrent control
#' for which the bias should be plotted.
#' @param out.sd.EC True sd of outcome in external control.
#' @param cov.CT List of covariates in concurrent treatment. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.CC List of covariates in concurrent control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.EC List of covariates in historical control. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cormat Matrix of correlation coefficients for outcome and covariates,
#' specified as Gaussian copula parameter.
#' @param method List of information borrowing method. \code{"noborrow"} uses
#' the concurrent data only. \code{"fullborrow"} uses the external control data
#' without discounting. \code{"cauchy"} uses the commensurate prior to
#' dynamically borrow the external control data, and the commensurate parameter
#' is assumed to follow half cauchy distribution. \code{"normal"} uses the
#' commensurate prior to dynamically borrow the external control data, and the
#' commensurate parameter is assumed to follow half normal distribution.
#' \code{"cauchy"} and \code{"normal"} require to specify the scale parameter of
#' half cauchy and half normal distribution respectively.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=2}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=4000}.
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.025}.
#' @param nsim Number of simulated trials. The default value is \code{nsim=10}.
#' @return
#' \item{reject}{\code{TRUE} when significant; otherwise \code{FALSE}.}
#' \item{theta}{Posterior mean, median, and sd of mean difference.}
#' @examples
#' n.CT  <- 100
#' n.CC  <- 100
#' n.EC  <- 100
#'
#' out.mean.CT <- 0
#' out.sd.CT   <- 1
#' out.mean.CC <- 0
#' out.sd.CC   <- 1
#' driftdiff   <- c(-0.2,0.0,0.2)
#' out.sd.EC   <- 1
#'
#' cov.CT <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cov.CC <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cov.EC <- list(list(dist="norm", mean=0,sd=1),
#'                list(dist="binom",prob=0.4))
#'
#' cormat <- rbind(c(  1,0.1,0.1),
#'                 c(0.1,  1,0.1),
#'                 c(0.1,0.1,  1))
#'
#' method <- list(list(prior="noborrow"),
#'                list(prior="fullborrow"),
#'                list(prior="cauchy",scale=0.5),
#'                list(prior="normal",scale=0.5))
#'
#' commensurate.cont(
#'   n.CT=n.CT, n.CC=n.CC, n.EC=n.EC,
#'   out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
#'   out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
#'   driftdiff=driftdiff, out.sd.EC=out.sd.EC,
#'   cov.CT=cov.CT, cov.CC=cov.CC, cov.EC=cov.EC, cormat=cormat,
#'   method=method)
#' @import rstan
#' @export

commensurate.cont <- function(
  n.CT, n.CC, n.EC,
  out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
  cov.CT, cov.CC, cov.EC, cormat, method,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim=10)
{
  ncov  <- length(cov.CT)
  out.mean.EC <- driftdiff+out.mean.CC
  nmethod     <- length(method)

  reject <- array(0,dim=c(nsim,nmethod))
  theta  <- array(0,dim=c(nsim,nmethod,3))

  for(ss in 1:nsim){

    marg.CT <- list(list(dist="norm",parm=list(mean=out.mean.CT,sd=out.sd.CT)))
    marg.CC <- list(list(dist="norm",parm=list(mean=out.mean.CC,sd=out.sd.CC)))
    marg.EC <- list(list(dist="norm",parm=list(mean=out.mean.EC,sd=out.sd.EC)))

    for(i in 1:ncov){
      if(cov.CT[[i]]$dist=="norm"){
        marg.CT <- append(marg.CT,list(list(dist=cov.CT[[i]]$dist,parm=list(mean=cov.CT[[i]]$mean,sd=cov.CT[[i]]$sd))))
        marg.CC <- append(marg.CC,list(list(dist=cov.CC[[i]]$dist,parm=list(mean=cov.CC[[i]]$mean,sd=cov.CC[[i]]$sd))))
        marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(mean=cov.EC[[i]]$mean,sd=cov.EC[[i]]$sd))))
      }else if(cov.CT[[i]]$dist=="binom"){
        marg.CT <- append(marg.CT,list(list(dist=cov.CT[[i]]$dist,parm=list(size=1,prob=cov.CT[[i]]$prob))))
        marg.CC <- append(marg.CC,list(list(dist=cov.CC[[i]]$dist,parm=list(size=1,prob=cov.CC[[i]]$prob))))
        marg.EC <- append(marg.EC,list(list(dist=cov.EC[[i]]$dist,parm=list(size=1,prob=cov.EC[[i]]$prob))))
      }
    }

    cvec <- cormat[lower.tri(cormat)]
    data.CT <- datagen(margdist=marg.CT,corvec=cvec,nsim=n.CT)
    data.CC <- datagen(margdist=marg.CC,corvec=cvec,nsim=n.CC)
    data.EC <- datagen(margdist=marg.EC,corvec=cvec,nsim=n.EC)

    for(i in 1:nmethod){

      if(method[[i]]$prior=="noborrow"){

        dat <- list(
          nCT = n.CT,
          nCC = n.CC,
          p   = ncov,
          yCT = data.CT[,1],
          yCC = data.CC[,1],
          xCT = data.CT[,-1],
          xCC = data.CC[,-1])

        mcmc <- rstan::sampling(stanmodels$ContNoborrow,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0)
        mcmc.sample <- rstan::extract(mcmc)

      }else if(method[[i]]$prior=="fullborrow"){

        dat <- list(
          nCT = n.CT,
          nCC = n.CC,
          nEC = n.EC,
          p   = ncov,
          yCT = data.CT[,1],
          yCC = data.CC[,1],
          yEC = data.EC[,1],
          xCT = data.CT[,-1],
          xCC = data.CC[,-1],
          xEC = data.EC[,-1])

        mcmc <- rstan::sampling(stanmodels$ContFullborrow,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0)
        mcmc.sample <- rstan::extract(mcmc)

      }else if(method[[i]]$prior=="cauchy"){

        dat <- list(
          nCT = n.CT,
          nCC = n.CC,
          nEC = n.EC,
          p   = ncov,
          yCT = data.CT[,1],
          yCC = data.CC[,1],
          yEC = data.EC[,1],
          xCT = data.CT[,-1],
          xCC = data.CC[,-1],
          xEC = data.EC[,-1],
          scale = method[[i]]$scale)

        mcmc <- rstan::sampling(stanmodels$ContCauchy,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0)
        mcmc.sample <- rstan::extract(mcmc)

      }else if(method[[i]]$prior=="normal"){

        dat <- list(
          nCT = n.CT,
          nCC = n.CC,
          nEC = n.EC,
          p   = ncov,
          yCT = data.CT[,1],
          yCC = data.CC[,1],
          yEC = data.EC[,1],
          xCT = data.CT[,-1],
          xCC = data.CC[,-1],
          xEC = data.EC[,-1],
          scale = method[[i]]$scale)

        mcmc <- rstan::sampling(stanmodels$ContNormal,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0)
        mcmc.sample <- rstan::extract(mcmc)

      }

      if(alternative=="greater"){
        postprob <- mean(mcmc.sample$theta>0)
      }else if(alternative=="less"){
        postprob <- mean(mcmc.sample$theta<0)
      }
      reject[ss,i] <- (postprob>(1-sig.level))
      theta[ss,i,] <- c(mean(mcmc.sample$theta),median(mcmc.sample$theta),sd(mcmc.sample$theta))
    }}

  return(list(reject=reject,theta=theta))
}

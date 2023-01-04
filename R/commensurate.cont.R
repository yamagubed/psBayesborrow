
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
#'   cov.C, cov.cor.C, cov.effect.C,
#'   cov.E, cov.cor.E, cov.effect.E,
#'   method, chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025, nsim)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param n.EC Number of patients in external control.
#' @param out.mean.CT True mean of outcome in concurrent treatment.
#' @param out.sd.CT True sd of outcome in concurrent treatment.
#' @param out.mean.CC True mean of outcome in concurrent control.
#' @param out.sd.CC True as of outcome in concurrent control.
#' @param driftdiff Difference between external control and concurrent control
#' for which the bias should be plotted.
#' @param cov.C List of covariates in concurrent data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.C Matrix of correlation coefficients for covariates in
#' concurrent data, specified as Gaussian copula parameter.
#' @param cov.effect.C Vector of covariate effects in concurrent data,
#' specified as mean difference.
#' @param cov.E List of covariates in external data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.E Matrix of correlation coefficients for covariates in
#' external data, specified as Gaussian copula parameter.
#' @param cov.effect.E Vector of covariate effects in external data,
#' specified as mean difference.
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
#' n.CT  <- 50
#' n.CC  <- 50
#' n.EC  <- 50
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
#' cov.E <- list(list(dist="norm",mean=0,sd=1),
#'               list(dist="binom",prob=0.4))
#'
#' cov.cor.E <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.E <- c(0.1,0.1)
#'
#' method <- list(list(prior="noborrow"),
#'                list(prior="normal",scale=0.5))
#'
#' nsim <- 5
#'
#' commensurate.cont(
#'   n.CT=n.CT, n.CC=n.CC, n.EC=n.EC,
#'   out.mean.CT=out.mean.CT, out.sd.CT=out.sd.CT,
#'   out.mean.CC=out.mean.CC, out.sd.CC=out.sd.CC,
#'   driftdiff=driftdiff, out.sd.EC=out.sd.EC,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.E=cov.E, cov.cor.E=cov.cor.E, cov.effect.E=cov.effect.E,
#'   method=method, nsim=nsim)
#' @import rstan
#' @export

commensurate.cont <- function(
  n.CT, n.CC, n.EC,
  out.mean.CT, out.sd.CT, out.mean.CC, out.sd.CC, driftdiff, out.sd.EC,
  cov.C, cov.cor.C, cov.effect.C,
  cov.E, cov.cor.E, cov.effect.E,
  method, chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim)
{
  ncov  <- length(cov.C)
  out.mean.EC <- driftdiff+out.mean.CC
  nmethod     <- length(method)

  marg.C <- NULL
  marg.E <- NULL
  mean.C <- NULL
  mean.E <- NULL

  for(i in 1:ncov){
    if(cov.C[[i]]$dist=="norm"){
      marg.C <- append(marg.C,list(list(dist=cov.C[[i]]$dist,parm=list(mean=cov.C[[i]]$mean,sd=cov.C[[i]]$sd))))
      marg.E <- append(marg.E,list(list(dist=cov.E[[i]]$dist,parm=list(mean=cov.E[[i]]$mean,sd=cov.E[[i]]$sd))))

      mean.C <- c(mean.C,cov.C[[i]]$mean)
      mean.E <- c(mean.E,cov.E[[i]]$mean)
    }else if(cov.C[[i]]$dist=="binom"){
      marg.C <- append(marg.C,list(list(dist=cov.C[[i]]$dist,parm=list(size=1,prob=cov.C[[i]]$prob))))
      marg.E <- append(marg.E,list(list(dist=cov.E[[i]]$dist,parm=list(size=1,prob=cov.E[[i]]$prob))))

      mean.C <- c(mean.C,cov.C[[i]]$prob)
      mean.E <- c(mean.E,cov.E[[i]]$prob)
    }
  }

  int.C   <- out.mean.CC-sum(mean.C*cov.effect.C)
  int.E   <- out.mean.EC-sum(mean.E*cov.effect.E)
  t.theta <- out.mean.CT-out.mean.CC

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.E  <- cov.cor.E[lower.tri(cov.cor.E)]

  reject <- array(0,dim=c(nsim,nmethod))
  theta  <- array(0,dim=c(nsim,nmethod,3))

  for(ss in 1:nsim){

    data.cov.CT <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CT)
    data.cov.CC <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CC)
    data.cov.EC <- datagen(margdist=marg.E,corvec=cvec.E,nsim=n.EC)

    mu.CT <- int.C+t.theta+apply(data.cov.CT,1,function(x){sum(x*cov.effect.C)})
    mu.CC <- int.C        +apply(data.cov.CC,1,function(x){sum(x*cov.effect.C)})
    mu.EC <- int.E        +apply(data.cov.EC,1,function(x){sum(x*cov.effect.E)})

    data.CT <- cbind(rnorm(n.CT,mean=mu.CT,sd=out.sd.CT),data.cov.CT)
    data.CC <- cbind(rnorm(n.CC,mean=mu.CC,sd=out.sd.CC),data.cov.CC)
    data.EC <- cbind(rnorm(n.EC,mean=mu.EC,sd=out.sd.EC),data.cov.EC)

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

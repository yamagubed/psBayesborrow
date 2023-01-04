
#' Simulation study of hybrid control design with commensurate power prior:
#' time-to-event outcome
#'
#' Bayesian dynamic information borrowing with commensurate power prior is
#' implemented for hybrid control design, where the concurrent control is
#' augmented by historical control. No borrowing and full borrowing are also
#' implemented. The time-to-event outcome is applicable.
#' @usage
#' commensurate.t2e <- function(
#'   n.CT, n.CC, nevent.C, n.EC, nevent.E, accrual,
#'   out.mevent.CT, out.mevent.CC, driftHR,
#'   cov.C, cov.cor.C, cov.effect.C,
#'   cov.E, cov.cor.E, cov.effect.E,
#'   method, chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025, nsim)
#' @param n.CT Number of patients in concurrent treatment.
#' @param n.CC Number of patients in concurrent control.
#' @param nevent.C Number of events in concurrent data.
#' @param n.EC Number of patients in external control.
#' @param nevent.E Number of events in external control.
#' @param accrual Accrual rate (number of enrolled patients per month).
#' @param out.mevent.CT True median time to event in concurrent treatment.
#' @param out.mevent.CC True median time to event in concurrent control.
#' @param driftHR Hazard ratio between external control and concurrent control
#' for which the bias should be plotted.
#' @param cov.C List of covariates in concurrent data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.C Matrix of correlation coefficients for covariates in
#' concurrent data, specified as Gaussian copula parameter.
#' @param cov.effect.C Vector of covariate effects in concurrent data,
#' specified as hazard ratio.
#' @param cov.E List of covariates in external data. Distribution and
#' and its parameters need to be specified for each covariate.
#' @param cov.cor.E Matrix of correlation coefficients for covariates in
#' external data, specified as Gaussian copula parameter.
#' @param cov.effect.E Vector of covariate effects in external data,
#' specified as hazard ratio.
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
#' \item{theta}{Posterior mean, median, and sd of log hazard ratio.}
#' @examples
#' n.CT      <- 50
#' n.CC      <- 50
#' nevent.C  <- 80
#' n.EC      <- 50
#' nevent.E  <- 40
#' accrual   <- 16
#'
#' out.mevent.CT <- 6
#' out.mevent.CC <- 6
#' driftHR       <- 1
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
#' commensurate.t2e(
#'   n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'   n.EC=n.EC, nevent.E=nevent.E, accrual=accrual,
#'   out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.E=cov.E, cov.cor.E=cov.cor.E, cov.effect.E=cov.effect.E,
#'   method=method, alternative="less", nsim=nsim)
#' @import rstan
#' @export

commensurate.t2e <- function(
  n.CT, n.CC, nevent.C, n.EC, nevent.E, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C, cov.cor.C, cov.effect.C,
  cov.E, cov.cor.E, cov.effect.E,
  method, chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim)
{
  ncov          <- length(cov.C)
  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  out.lambda.EC <- out.lambda.CC*driftHR
  nmethod       <- length(method)

  lcov.effect.C <- (-log(cov.effect.C))
  lcov.effect.E <- (-log(cov.effect.E))

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

  int.C   <- log(1/out.lambda.CC)-sum(mean.C*lcov.effect.C)
  int.E   <- log(1/out.lambda.EC)-sum(mean.E*lcov.effect.E)
  t.theta <- log(out.lambda.CC/out.lambda.CT)

  cvec.C  <- cov.cor.C[lower.tri(cov.cor.C)]
  cvec.E  <- cov.cor.E[lower.tri(cov.cor.E)]

  reject <- array(0,dim=c(nsim,nmethod))
  theta  <- array(0,dim=c(nsim,nmethod,3))

  for(ss in 1:nsim){

    data.cov.CT <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CT)
    data.cov.CC <- datagen(margdist=marg.C,corvec=cvec.C,nsim=n.CC)
    data.cov.EC <- datagen(margdist=marg.E,corvec=cvec.E,nsim=n.EC)

    sigma.CT <- exp(int.C+t.theta+apply(data.cov.CT,1,function(x){sum(x*lcov.effect.C)}))
    sigma.CC <- exp(int.C        +apply(data.cov.CC,1,function(x){sum(x*lcov.effect.C)}))
    sigma.EC <- exp(int.E        +apply(data.cov.EC,1,function(x){sum(x*lcov.effect.E)}))

    data.CT <- cbind(rweibull(n.CT,shape=1,scale=sigma.CT),data.cov.CT)
    data.CC <- cbind(rweibull(n.CC,shape=1,scale=sigma.CC),data.cov.CC)
    data.EC <- cbind(rweibull(n.EC,shape=1,scale=sigma.EC),data.cov.EC)

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

    enroll.period.E <- floor(n.EC/accrual)

    enroll.time.EC <- NULL
    for(i in 1:enroll.period.E){
      e.st <- i-1
      e.en <- i
      enroll.time.EC <- c(enroll.time.EC,runif(accrual,e.st,e.en))
    }
    enroll.time.EC <- c(enroll.time.EC,runif(n.EC-length(enroll.time.EC),enroll.period.E,enroll.period.E+1))

    obs.time.CT <- data.CT[,1]+enroll.time.CT
    obs.time.CC <- data.CC[,1]+enroll.time.CC
    obs.time.EC <- data.EC[,1]+enroll.time.EC

    last.sub.C <- sort(c(obs.time.CT,obs.time.CC))[nevent.C]
    last.sub.E <- sort(obs.time.EC)[nevent.E]

    censor.CT <- as.numeric(obs.time.CT>last.sub.C)
    censor.CC <- as.numeric(obs.time.CC>last.sub.C)
    censor.EC <- as.numeric(obs.time.EC>last.sub.E)

    cenf <- (obs.time.CT>last.sub.C)
    data.CT[cenf,1] <- data.CT[cenf,1]-(obs.time.CT[cenf]-last.sub.C)

    cenf <- (obs.time.CC>last.sub.C)
    data.CC[cenf,1] <- data.CC[cenf,1]-(obs.time.CC[cenf]-last.sub.C)

    cenf <- (obs.time.EC>last.sub.E)
    data.EC[cenf,1] <- data.EC[cenf,1]-(obs.time.EC[cenf]-last.sub.E)

    censor.CT <- censor.CT[data.CT[,1]>0]
    data.CT   <- data.CT[data.CT[,1]>0,]

    censor.CC <- censor.CC[data.CC[,1]>0]
    data.CC   <- data.CC[data.CC[,1]>0,]

    censor.EC <- censor.EC[data.EC[,1]>0]
    data.EC   <- data.EC[data.EC[,1]>0,]

    for(i in 1:nmethod){

      if(sum(censor.CC==1)>0){

        if(method[[i]]$prior=="noborrow"){

          dat <- list(
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nCC_c = sum(censor.CC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yCC_c = as.array(data.CC[censor.CC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xCC_c = matrix(data.CC[censor.CC==1,-1],sum(censor.CC==1),ncov))

          mcmc <- rstan::sampling(stanmodels$T2ENoborrow,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nCC_c = sum(censor.CC==1),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yCC_c = as.array(data.CC[censor.CC==1,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xCC_c = matrix(data.CC[censor.CC==1,-1],sum(censor.CC==1),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov))

          mcmc <- rstan::sampling(stanmodels$T2EFullborrow,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nCC_c = sum(censor.CC==1),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yCC_c = as.array(data.CC[censor.CC==1,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xCC_c = matrix(data.CC[censor.CC==1,-1],sum(censor.CC==1),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
            scale = method[[i]]$scale)

          mcmc <- rstan::sampling(stanmodels$T2ECauchy,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nCC_c = sum(censor.CC==1),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yCC_c = as.array(data.CC[censor.CC==1,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xCC_c = matrix(data.CC[censor.CC==1,-1],sum(censor.CC==1),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
            scale = method[[i]]$scale)

          mcmc <- rstan::sampling(stanmodels$T2ENormal,
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

      }else if(sum(censor.CC==1)==0){

        if(method[[i]]$prior=="noborrow"){

          dat <- list(
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov))

          mcmc <- rstan::sampling(stanmodels$T2ENoborrowC0,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov))

          mcmc <- rstan::sampling(stanmodels$T2EFullborrowC0,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
            scale = method[[i]]$scale)

          mcmc <- rstan::sampling(stanmodels$T2ECauchyC0,
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
            nCT_o = sum(censor.CT==0),
            nCT_c = sum(censor.CT==1),
            nCC_o = sum(censor.CC==0),
            nEC_o = sum(censor.EC==0),
            nEC_c = sum(censor.EC==1),
            p     = ncov,
            yCT_o = as.array(data.CT[censor.CT==0,1]),
            yCT_c = as.array(data.CT[censor.CT==1,1]),
            yCC_o = as.array(data.CC[censor.CC==0,1]),
            yEC_o = as.array(data.EC[censor.EC==0,1]),
            yEC_c = as.array(data.EC[censor.EC==1,1]),
            xCT_o = matrix(data.CT[censor.CT==0,-1],sum(censor.CT==0),ncov),
            xCT_c = matrix(data.CT[censor.CT==1,-1],sum(censor.CT==1),ncov),
            xCC_o = matrix(data.CC[censor.CC==0,-1],sum(censor.CC==0),ncov),
            xEC_o = matrix(data.EC[censor.EC==0,-1],sum(censor.EC==0),ncov),
            xEC_c = matrix(data.EC[censor.EC==1,-1],sum(censor.EC==1),ncov),
            scale = method[[i]]$scale)

          mcmc <- rstan::sampling(stanmodels$T2ENormalC0,
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

      }

      loghr <- (-mcmc.sample$alpha*mcmc.sample$theta)

      if(alternative=="greater"){
        postprob <- mean(loghr>0)
      }else if(alternative=="less"){
        postprob <- mean(loghr<0)
      }
      reject[ss,i] <- (postprob>(1-sig.level))
      theta[ss,i,] <- c(mean(loghr),median(loghr),sd(loghr))
    }}

  return(list(reject=reject,theta=theta))
}

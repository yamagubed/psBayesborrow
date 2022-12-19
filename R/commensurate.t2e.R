
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
#'   cov.CT, cov.CC, cov.EC, cormat, method=method,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025, nsim=10)
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
#' \item{theta}{Posterior mean, median, and sd of log hazard ratio.}
#' @examples
#' n.CT      <- 60
#' n.CC      <- 60
#' nevent.C  <- 100
#' n.EC      <- 120
#' nevent.E  <- 100
#' accrual   <- 16
#'
#' out.mevent.CT <- 6
#' out.mevent.CC <- 6
#' driftHR       <- 1.0
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
#' commensurate.t2e(
#'   n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'   n.EC=n.EC, nevent.E=nevent.E, accrual=accrual,
#'   out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
#'   cov.CT=cov.CT, cov.CC=cov.CC, cov.EC=cov.EC, cormat=cormat,
#'   method=method)
#' @import rstan
#' @export

commensurate.t2e <- function(
  n.CT, n.CC, nevent.C, n.EC, nevent.E, accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.CT, cov.CC, cov.EC, cormat, method=method,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim=10)
{
  ncov          <- length(cov.CT)
  out.lambda.CT <- log(2)/out.mevent.CT
  out.lambda.CC <- log(2)/out.mevent.CC
  out.lambda.EC <- out.lambda.CC/driftHR
  nmethod       <- length(method)

  reject <- array(0,dim=c(nsim,nmethod))
  theta  <- array(0,dim=c(nsim,nmethod,3))

  for(ss in 1:nsim){

    marg.CT <- list(list(dist="exp",parm=list(rate=out.lambda.CT)))
    marg.CC <- list(list(dist="exp",parm=list(rate=out.lambda.CC)))
    marg.EC <- list(list(dist="exp",parm=list(rate=out.lambda.EC)))

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

    enroll.period.C <- floor((n.CT+n.CC)/accrual)

    enroll.time.CT <- NULL
    enroll.time.CC <- NULL
    for(i in 1:enroll.period.C){
      tmp <- runif(1)
      if(tmp<0.5){
        tmp.n.CT <- floor(accrual*(n.CT/(n.CT+n.CC)))
      }else{
        tmp.n.CT <- ceiling(accrual*(n.CT/(n.CT+n.CC)))
      }
      tmp.n.CC <- accrual-tmp.n.CT

      e.st <- i-1
      e.en <- i
      enroll.time.CT <- c(enroll.time.CT,runif(tmp.n.CT,e.st,e.en))
      enroll.time.CC <- c(enroll.time.CC,runif(tmp.n.CC,e.st,e.en))
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

      if(method[[i]]$prior=="noborrow"){

        dat <- list(
          nCT_o = sum(censor.CT==0),
          nCT_c = sum(censor.CT==1),
          nCC_o = sum(censor.CC==0),
          nCC_c = sum(censor.CC==1),
          p     = ncov,
          yCT_o = data.CT[censor.CT==0,1],
          yCT_c = data.CT[censor.CT==1,1],
          yCC_o = data.CC[censor.CC==0,1],
          yCC_c = data.CC[censor.CC==1,1],
          xCT_o = data.CT[censor.CT==0,-1],
          xCT_c = data.CT[censor.CT==1,-1],
          xCC_o = data.CC[censor.CC==0,-1],
          xCC_c = data.CC[censor.CC==1,-1])

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
          yCT_o = data.CT[censor.CT==0,1],
          yCT_c = data.CT[censor.CT==1,1],
          yCC_o = data.CC[censor.CC==0,1],
          yCC_c = data.CC[censor.CC==1,1],
          yEC_o = data.EC[censor.EC==0,1],
          yEC_c = data.EC[censor.EC==1,1],
          xCT_o = data.CT[censor.CT==0,-1],
          xCT_c = data.CT[censor.CT==1,-1],
          xCC_o = data.CC[censor.CC==0,-1],
          xCC_c = data.CC[censor.CC==1,-1],
          xEC_o = data.EC[censor.EC==0,-1],
          xEC_c = data.EC[censor.EC==1,-1])

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
          yCT_o = data.CT[censor.CT==0,1],
          yCT_c = data.CT[censor.CT==1,1],
          yCC_o = data.CC[censor.CC==0,1],
          yCC_c = data.CC[censor.CC==1,1],
          yEC_o = data.EC[censor.EC==0,1],
          yEC_c = data.EC[censor.EC==1,1],
          xCT_o = data.CT[censor.CT==0,-1],
          xCT_c = data.CT[censor.CT==1,-1],
          xCC_o = data.CC[censor.CC==0,-1],
          xCC_c = data.CC[censor.CC==1,-1],
          xEC_o = data.EC[censor.EC==0,-1],
          xEC_c = data.EC[censor.EC==1,-1],
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
          yCT_o = data.CT[censor.CT==0,1],
          yCT_c = data.CT[censor.CT==1,1],
          yCC_o = data.CC[censor.CC==0,1],
          yCC_c = data.CC[censor.CC==1,1],
          yEC_o = data.EC[censor.EC==0,1],
          yEC_c = data.EC[censor.EC==1,1],
          xCT_o = data.CT[censor.CT==0,-1],
          xCT_c = data.CT[censor.CT==1,-1],
          xCC_o = data.CC[censor.CC==0,-1],
          xCC_c = data.CC[censor.CC==1,-1],
          xEC_o = data.EC[censor.EC==0,-1],
          xEC_c = data.EC[censor.EC==1,-1],
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

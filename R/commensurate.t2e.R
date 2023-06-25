
#' Commensurate power prior: time-to-event outcome
#'
#' Commensurate power prior is implemented. No borrowing and full borrowing are
#' also implemented. The time-to-event outcome is applicable.
#' @usage
#' commensurate.t2e(
#'   indata, subjid.EC, method.borrow,
#'   chains=2, iter=4000, warmup=floor(iter/2), thin=1,
#'   alternative="greater", sig.level=0.025)
#' @param indata Dataset of a simulated trial, which is a data frame returned
#' from \code{trial.simulation.t2e}, \code{trial.simulation.bin}, or
#' \code{trial.simulation.cont}.
#' @param subjid.EC Subject ID of external control.
#' @param method.borrow List of information borrowing method. \code{"noborrow"}
#' uses the concurrent data only. \code{"fullborrow"} uses the external control
#' data without discounting. \code{"cauchy"} uses the commensurate prior to
#' dynamically borrow the external control data, and the commensurate parameter
#' is assumed to follow half-Cauchy distribution. \code{"normal"} uses the
#' commensurate prior to dynamically borrow the external control data, and the
#' commensurate parameter is assumed to follow half-normal distribution.
#' \code{"cauchy"} and \code{"normal"} require to specify the scale parameter of
#' half-Cauchy and half-normal distribution respectively.
#' @param chains Number of Markov chains in MCMC sampling. The default value is
#' \code{chains=2}.
#' @param iter Number of iterations for each chain (including warmup) in MCMC
#' sampling. The default value is \code{iter=4000}.
#' @param warmup Number of warmup (burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.025}.
#' @return
#' \item{reject}{\code{TRUE} when significant; otherwise \code{FALSE}.}
#' \item{theta}{Posterior mean, median, and sd of log hazard ratio.}
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
#' indata <- trial.simulation.t2e(
#'   n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
#'   n.ECp=n.ECp, nevent.ECp=nevent.ECp, accrual=accrual,
#'   out.mevent.CT, out.mevent.CC, driftHR,
#'   cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
#'   cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC)
#'
#' n.EC <- 50
#'
#' method.whomatch <- "conc.treat"
#' method.matching <- "optimal"
#' method.psorder  <- NULL
#'
#' out.psmatch <- psmatch(
#'   indata=indata, n.EC=n.EC,
#'   method.whomatch=method.whomatch, method.matching=method.matching,
#'   method.psorder=method.psorder)
#'
#' subjid.EC <- out.psmatch$subjid.EC
#'
#' method.borrow <- list(list(prior="noborrow"),
#'                       list(prior="normal",scale=0.5))
#'
#' commensurate.t2e(
#'   indata=indata, subjid.EC=subjid.EC, method.borrow=method.borrow)
#' @import rstan
#' @export

commensurate.t2e <- function(
  indata, subjid.EC, method.borrow,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025)
{
  ncov      <- indata$ncov
  cov.lab   <- paste("X",1:ncov,sep="")

  data.outp <- indata$data
  data.out  <- rbind(data.outp[data.outp$study==1,],data.outp[subjid.EC,])

  censor.CT <- data.out[(data.out$study==1)&(data.out$treat==1),"censor"]
  data.CT   <- data.out[(data.out$study==1)&(data.out$treat==1),c("y",cov.lab)]

  censor.CC <- data.out[(data.out$study==1)&(data.out$treat==0),"censor"]
  data.CC   <- data.out[(data.out$study==1)&(data.out$treat==0),c("y",cov.lab)]

  censor.EC <- data.out[(data.out$study==0)&(data.out$treat==0),"censor"]
  data.EC   <- data.out[(data.out$study==0)&(data.out$treat==0),c("y",cov.lab)]

  nmethod    <- length(method.borrow)
  method.lab <- character(nmethod)

  reject <- data.frame(X0="reject")
  theta  <- data.frame(X0=c("mean","median","sd","lcri","ucri"))

  for(i in 1:nmethod){

    if(sum(censor.CC==1)>0){

      if(method.borrow[[i]]$prior=="noborrow"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xCC_c = matrix(as.matrix(data.CC[censor.CC==1,-1]),sum(censor.CC==1),ncov))

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

      }else if(method.borrow[[i]]$prior=="fullborrow"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xCC_c = matrix(as.matrix(data.CC[censor.CC==1,-1]),sum(censor.CC==1),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov))

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

      }else if(method.borrow[[i]]$prior=="cauchy"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xCC_c = matrix(as.matrix(data.CC[censor.CC==1,-1]),sum(censor.CC==1),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov),
          scale = method.borrow[[i]]$scale)

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

      }else if(method.borrow[[i]]$prior=="normal"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xCC_c = matrix(as.matrix(data.CC[censor.CC==1,-1]),sum(censor.CC==1),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov),
          scale = method.borrow[[i]]$scale)

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

      if(method.borrow[[i]]$prior=="noborrow"){

        dat <- list(
          nCT_o = sum(censor.CT==0),
          nCT_c = sum(censor.CT==1),
          nCC_o = sum(censor.CC==0),
          p     = ncov,
          yCT_o = as.array(data.CT[censor.CT==0,1]),
          yCT_c = as.array(data.CT[censor.CT==1,1]),
          yCC_o = as.array(data.CC[censor.CC==0,1]),
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov))

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

      }else if(method.borrow[[i]]$prior=="fullborrow"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov))

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

      }else if(method.borrow[[i]]$prior=="cauchy"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov),
          scale = method.borrow[[i]]$scale)

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

      }else if(method.borrow[[i]]$prior=="normal"){

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
          xCT_o = matrix(as.matrix(data.CT[censor.CT==0,-1]),sum(censor.CT==0),ncov),
          xCT_c = matrix(as.matrix(data.CT[censor.CT==1,-1]),sum(censor.CT==1),ncov),
          xCC_o = matrix(as.matrix(data.CC[censor.CC==0,-1]),sum(censor.CC==0),ncov),
          xEC_o = matrix(as.matrix(data.EC[censor.EC==0,-1]),sum(censor.EC==0),ncov),
          xEC_c = matrix(as.matrix(data.EC[censor.EC==1,-1]),sum(censor.EC==1),ncov),
          scale = method.borrow[[i]]$scale)

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

    loghr <- mcmc.sample$theta

    if(alternative=="greater"){
      postprob <- mean(loghr>0)
    }else if(alternative=="less"){
      postprob <- mean(loghr<0)
    }
    cri <- quantile(loghr,c(sig.level,1-sig.level))

    reject.v <- (postprob>(1-sig.level))
    reject   <- data.frame(reject,X1=reject.v)

    theta.v <- c(mean(loghr),median(loghr),sd(loghr),cri[[1]],cri[[2]])
    theta   <- data.frame(theta,X1=theta.v)

    mname <- method.borrow[[i]]$prior
    if((mname=="noborrow")|(mname=="fullborrow")){
      method.lab[i] <- mname
    }else if((mname=="cauchy")|(mname=="normal")){
      method.lab[i] <- paste(mname,method.borrow[[i]]$scale,sep="")
    }
  }

  colnames(reject) <- c("measure",method.lab)
  colnames(theta)  <- c("measure",method.lab)

  return(list(reject=reject,theta=theta,method.lab=method.lab))
}

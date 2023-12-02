
#' Commensurate power prior: binary outcome
#'
#' Commensurate power prior is implemented. No borrowing and full borrowing are
#' also implemented. The binary outcome is applicable.
#' @usage
#' commensurate.bin(
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
#' @param warmup Number of warmup (aka burnin) iterations per chain in MCMC
#' sampling. The default value is \code{warmup=floor(iter/2)}.
#' @param thin Period for saving samples in MCMC sampling. The default value
#' is \code{thin=1}.
#' @param alternative Alternative hypothesis to be tested ("greater" or "less").
#' The default value is \code{alternative="greater"}.
#' @param sig.level Significance level. The default value is
#' \code{sig.level=0.025}.
#' @return
#' \item{reject}{\code{TRUE} when significant; otherwise \code{FALSE}.}
#' \item{theta}{Posterior mean, median, and sd of log odds ratio.}
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
#' indata <- trial.simulation.bin(
#'   n.CT=n.CT, n.CC=n.CC, n.ECp=n.ECp,
#'   out.prob.CT=out.prob.CT, out.prob.CC=out.prob.CC, driftOR=driftOR,
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
#' commensurate.bin(
#'   indata=indata, subjid.EC=subjid.EC, method.borrow=method.borrow)
#' @import rstan
#' @export

commensurate.bin <- function(
  formula, data, method.borrow,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025,
  seed=sample.int(.Machine$integer.max,1))
{
  mf      <- model.frame(formula=formula,data=data)
  cov.lab <- attr(attr(mf,"terms"),"term.labels")
  y       <- model.response(mf)

  if(length(data$study)==0){
    stop("Dataset must contain a variable of study.")

  }else if(length(data$treat)==0){
    stop("Dataset must contain a variable of treat.")

  }else{

    ncov     <- length(cov.lab)
    data.out <- data.frame(y      = as.vector(y),
                           study  = data$study,
                           treat  = data$treat,
                           data[,cov.lab,drop=FALSE])

    data.CT <- data.out[(data.out$study==1)&(data.out$treat==1),c("y",cov.lab)]
    data.CC <- data.out[(data.out$study==1)&(data.out$treat==0),c("y",cov.lab)]
    data.EC <- data.out[(data.out$study==0)&(data.out$treat==0),c("y",cov.lab)]

    nmethod    <- length(method.borrow)
    method.lab <- character(nmethod)

    reject <- data.frame(X0="reject")
    theta  <- data.frame(X0=c("mean","median","sd","lcri","ucri"))

    stan.obj <- NULL

    for(i in 1:nmethod){

      if(method.borrow[[i]]$prior=="noborrow"){

        dat <- list(
          nCT = nrow(data.CT),
          nCC = nrow(data.CC),
          p   = ncov,
          yCT = as.array(data.CT[,"y"]),
          yCC = as.array(data.CC[,"y"]),
          xCT = as.matrix(data.CT[,cov.lab,drop=FALSE]),
          xCC = as.matrix(data.CC[,cov.lab,drop=FALSE]))

        mcmc <- rstan::sampling(stanmodels$BinNoborrow,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0,
                                seed          = seed)
        mcmc.sample <- rstan::extract(mcmc)

    }else if(method.borrow[[i]]$prior=="fullborrow"){

        dat <- list(
          nCT = nrow(data.CT),
          nCC = nrow(data.CC),
          nEC = nrow(data.EC),
          p   = ncov,
          yCT = as.array(data.CT[,"y"]),
          yCC = as.array(data.CC[,"y"]),
          yEC = as.array(data.EC[,"y"]),
          xCT = as.matrix(data.CT[,cov.lab,drop=FALSE]),
          xCC = as.matrix(data.CC[,cov.lab,drop=FALSE]),
          xEC = as.matrix(data.EC[,cov.lab,drop=FALSE]))

        mcmc <- rstan::sampling(stanmodels$BinFullborrow,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0,
                                seed          = seed)
        mcmc.sample <- rstan::extract(mcmc)

      }else if(method.borrow[[i]]$prior=="cauchy"){

        dat <- list(
          nCT   = nrow(data.CT),
          nCC   = nrow(data.CC),
          nEC   = nrow(data.EC),
          p     = ncov,
          yCT   = as.array(data.CT[,"y"]),
          yCC   = as.array(data.CC[,"y"]),
          yEC   = as.array(data.EC[,"y"]),
          xCT   = as.matrix(data.CT[,cov.lab,drop=FALSE]),
          xCC   = as.matrix(data.CC[,cov.lab,drop=FALSE]),
          xEC   = as.matrix(data.EC[,cov.lab,drop=FALSE]),
          scale = method.borrow[[i]]$scale)

        mcmc <- rstan::sampling(stanmodels$BinCauchy,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0,
                                seed          = seed)
        mcmc.sample <- rstan::extract(mcmc)

      }else if(method.borrow[[i]]$prior=="normal"){

        dat <- list(
          nCT   = nrow(data.CT),
          nCC   = nrow(data.CC),
          nEC   = nrow(data.EC),
          p     = ncov,
          yCT   = as.array(data.CT[,"y"]),
          yCC   = as.array(data.CC[,"y"]),
          yEC   = as.array(data.EC[,"y"]),
          xCT   = as.matrix(data.CT[,cov.lab,drop=FALSE]),
          xCC   = as.matrix(data.CC[,cov.lab,drop=FALSE]),
          xEC   = as.matrix(data.EC[,cov.lab,drop=FALSE]),
          scale = method.borrow[[i]]$scale)

        mcmc <- rstan::sampling(stanmodels$BinNormal,
                                data          = dat,
                                chains        = chains,
                                iter          = iter,
                                warmup        = warmup,
                                thin          = thin,
                                show_messages = FALSE,
                                cores         = 1,
                                refresh       = 0,
                                seed          = seed)
        mcmc.sample <- rstan::extract(mcmc)

      }

      logor <- mcmc.sample$theta

      if(alternative=="greater"){
        postprob <- mean(logor>0)
      }else if(alternative=="less"){
        postprob <- mean(logor<0)
      }
      cri <- quantile(logor,c(sig.level,1-sig.level))

      reject.v <- (postprob>(1-sig.level))
      reject   <- data.frame(reject,X1=reject.v)

      theta.v <- c(mean(logor),median(logor),sd(logor),cri[[1]],cri[[2]])
      theta   <- data.frame(theta,X1=theta.v)

      mname <- method.borrow[[i]]$prior
      if((mname=="noborrow")|(mname=="fullborrow")){
        method.lab[i] <- mname
      }else if((mname=="cauchy")|(mname=="normal")){
        method.lab[i] <- paste(mname,method.borrow[[i]]$scale,sep="")
      }

      stan.obj <- append(stan.obj,list(mcmc))
    }

    colnames(reject) <- c("measure",method.lab)
    colnames(theta)  <- c("measure",method.lab)
    names(stan.obj)  <- method.lab

    return(list(reject=reject,theta=theta,stan.obj=stan.obj))
  }
}

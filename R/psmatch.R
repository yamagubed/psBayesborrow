
#' Propensity score matching
#'
#' Propensity score matching is implemented using methods specified.
#' @usage
#' psmatch(
#'   indata, n.EC,
#'   method.psest="glm", method.pslink="logit",
#'   method.whomatch, method.matching, method.psorder)
#' @param indata Dataset of a simulated trial, which is a data frame returned
#' from \code{trial.simulation.t2e}, \code{trial.simulation.bin}, or
#' \code{trial.simulation.cont}.
#' @param n.EC Number of patients in external control.
#' @param method.psest Method of estimating the propensity score. Allowable
#' options include, for example, \code{"glm"} for generalized linear model
#' (e.g., logistic regression); \code{"gam"} for generalized additive model;
#' \code{"gbm"} for generalized boosted model; \code{"lasso"} for lasso
#' regression; \code{"rpart"} for classification tree. The default value is
#' \code{method.psest="glm"}.
#' @param method.pslink Link function used in estimating the propensity score.
#' Allowable options depend on the specific \code{method.psest} value specified.
#' The default value is \code{method.pslink="logit"}, which, along with
#' \code{method.psest="glm"}, identifies the default method as logistic
#' regression.
#' @param method.whomatch Options of who to match. Allowable options include
#' \code{conc.contl} matching concurrent control to external control pool;
#' \code{conc.treat} matching concurrent treatment to external control pool;
#' \code{conc.all} matching concurrent treatment plus concurrent control to
#' external control pool; \code{treat2contl} matching concurrent treatment to
#' concurrent control plus external control pool.
#' @param method.matching Matching method. Allowable options include
#' \code{"optimal"} for optimal matching; \code{"nearest"} for nearest neighbor
#' matching without replacement.
#' @param method.psorder Order that the matching takes place. Allowable options
#' include \code{"largest"}, where matching takes place in descending order of
#' propensity score; \code{"smallest"}, where matching takes place in ascending
#' order of propensity score; \code{"random"}, where matching takes place in a
#' random order; \code{"data"}, where matching takes place based on the order
#' of units in the data.
#' @return
#' \item{subjid.EC}{Subject ID of external control.}
#' \item{data.cov.ps}{Dataset of a simulated trial with estimated propensity
#' score.}
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
#' psmatch(
#'   indata=indata, n.EC=n.EC,
#'   method.whomatch=method.whomatch, method.matching=method.matching,
#'   method.psorder=method.psorder)
#' @import MatchIt
#' @export

psmatch <- function(
  indata, n.EC,
  method.psest="glm", method.pslink="logit",
  method.whomatch, method.matching, method.psorder)
{
  cov.lab  <- paste("X",1:indata$ncov,sep="")
  data.cov <- data.frame(indata$data[,c("study","treat",cov.lab)],
                         conc.contl=0,conc.treat=0,conc.all=0,treat2contl=0)

  CT.flg <- (data.cov$study==1)&(data.cov$treat==1)
  CC.flg <- (data.cov$study==1)&(data.cov$treat==0)

  data.cov[CT.flg,"conc.contl"]  <- 99
  data.cov[CT.flg,"conc.treat"]  <- 1
  data.cov[CT.flg,"conc.all"]    <- 1
  data.cov[CT.flg,"treat2contl"] <- 1

  data.cov[CC.flg,"conc.contl"]  <- 1
  data.cov[CC.flg,"conc.treat"]  <- 99
  data.cov[CC.flg,"conc.all"]    <- 1
  data.cov[CC.flg,"treat2contl"] <- 0

  fit.ps <- MatchIt::matchit(study~.,
                             data     = data.cov[,c("study",cov.lab)],
                             method   = NULL,
                             distance = method.psest,
                             link     = method.pslink)
  data.cov.ps <- data.frame(data.cov,ps=as.vector(fit.ps$distance))

  match.pts  <- (data.cov.ps[,method.whomatch]!=99)
  data.match <- data.cov.ps[match.pts,c(method.whomatch,cov.lab)]
  colnames(data.match) <- c("match",cov.lab)

  if(method.whomatch=="conc.contl"){

    n.CC  <- sum(CC.flg)
    ratio <- n.EC/n.CC

    if(method.matching=="optimal"){

      fit.match <- MatchIt::matchit(match~.,
                                    data     = data.match,
                                    method   = method.matching,
                                    distance = data.cov.ps[match.pts,"ps"],
                                    ratio    = ratio)
      subjid.EC <- as.numeric(fit.match$match.matrix)

    }else if(method.matching=="nearest"){

      fit.match <- MatchIt::matchit(match~.,
                                    data     = data.match,
                                    method   = method.matching,
                                    distance = data.cov.ps[match.pts,"ps"],
                                    m.order  = method.psorder,
                                    ratio    = ratio)
      subjid.EC <- as.numeric(fit.match$match.matrix)
    }

  }else if((method.whomatch=="conc.treat")|(method.whomatch=="conc.all")){

    if(method.matching=="optimal"){

      fit.match <- MatchIt::matchit(match~.,
                                    data     = data.match,
                                    method   = method.matching,
                                    distance = data.cov.ps[match.pts,"ps"])
      subjid.EC <- sample(as.numeric(fit.match$match.matrix),n.EC)

    }else if(method.matching=="nearest"){

      ps.M1 <- data.cov.ps[data.cov.ps[,method.whomatch]==1,"ps"]
      ps.M0 <- data.cov.ps[data.cov.ps[,method.whomatch]==0,"ps"]
      ps.M  <- sort(abs(apply(as.matrix(ps.M1),1,function(x){x-ps.M0})))

      ncal <- n.EC
      while(1){
        fit.match <- MatchIt::matchit(match~.,
                                      data        = data.match,
                                      method      = method.matching,
                                      distance    = data.cov.ps[match.pts,"ps"],
                                      m.order     = method.psorder,
                                      caliper     = ps.M[ncal],
                                      std.caliper = FALSE)
        if(sum(!is.na(fit.match$match.matrix))>=n.EC){
          nearest <- fit.match$match.matrix
          break
        }else{
          ncal <- ncal+1
        }
      }
      subjid.EC <- sample(as.numeric(nearest[!is.na(nearest),]),n.EC)
    }

  }else if(method.whomatch=="treat2contl"){

    if(method.matching=="optimal"){

      fit.match <- MatchIt::matchit(match~.,
                                    data     = data.match,
                                    method   = method.matching,
                                    distance = data.cov.ps[match.pts,"ps"])
      subjid.C1 <- as.numeric(fit.match$match.matrix)
      subjid.C2 <- subjid.C1[data.cov.ps[subjid.C1,"study"]==0]
      if(length(subjid.C2)>=n.EC){
        subjid.EC <- sample(subjid.C2,n.EC)
      }else{
        subjid.EC <- subjid.C2
      }

    }else if(method.matching=="nearest"){

      ps.M1 <- data.cov.ps[data.cov.ps[,method.whomatch]==1,"ps"]
      ps.M0 <- data.cov.ps[data.cov.ps[,method.whomatch]==0,"ps"]
      ps.M  <- sort(abs(apply(as.matrix(ps.M1),1,function(x){x-ps.M0})))

      ncal <- n.EC
      while(1){
        fit.match <- MatchIt::matchit(match~.,
                                      data        = data.match,
                                      method      = method.matching,
                                      distance    = data.cov.ps[match.pts,"ps"],
                                      m.order     = method.psorder,
                                      caliper     = ps.M[ncal],
                                      std.caliper = FALSE)

        subjid.C1 <- fit.match$match.matrix
        subjid.C2 <- as.numeric(subjid.C1[!is.na(subjid.C1)])
        subjid.C3 <- subjid.C2[data.cov.ps[subjid.C2,"study"]==0]

        if(length(subjid.C3)>=n.EC){
          nearest <- subjid.C3
          break
        }else{
          ncal <- ncal+1
        }
      }
      subjid.EC <- sample(nearest,n.EC)
    }
  }

  return(list(subjid.EC=subjid.EC,data.cov.ps=data.cov.ps))
}

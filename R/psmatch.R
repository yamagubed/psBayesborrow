
#' Propensity score matching
#'
#' Propensity score matching is implemented using methods specified.
#' @usage
#' psmatch(
#'     formula, data, n.EC,
#'     method.psest="glm", method.pslink="logit",
#'     method.whomatch, method.matching, method.psorder)
#' @param formula formula.
#' @param data data.
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
#' cov.C <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'               list(dist="binom",prob=0.4,lab="cov2"))
#'
#' cov.cor.C <- rbind(c(  1,0.1),
#'                    c(0.1,  1))
#'
#' cov.effect.C <- c(0.1,0.1)
#'
#' cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
#'                list(dist="binom",prob=0.4,lab="cov2"))
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
#'   study~cov1+cov2, data=indata, n.EC=n.EC,
#'   method.whomatch=method.whomatch, method.matching=method.matching,
#'   method.psorder=method.psorder)
#' @import dplyr optmatch MatchIt e1071
#' @export

psmatch <- function(
  formula, data, n.EC,
  method.psest="glm", method.pslink="logit",
  method.whomatch, method.matching, method.psorder, n.boot=100)
{
  mf      <- model.frame(formula=formula,data=data)
  cov.lab <- attr(attr(mf,"terms"),"term.labels")

  if(length(data$study)==0){
    stop("Dataset must contain a variable of study.")

  }else if(length(data$treat)==0){
    stop("Dataset must contain a variable of treat.")

  }else{

    bootmatch.optimal <- function(nwps,dmat)
    {
      boot.data.match <- rbind(data.frame(match=1,ps=nwps),dmat[dmat$match==0,])
      bootmatching    <- MatchIt::matchit(match~.,
                                          data     = boot.data.match,
                                          method   = "optimal",
                                          distance = boot.data.match[,"ps"])
      subjid.EC  <- as.numeric(bootmatching$match.matrix)
      return(subjid.EC)
    }

    bootmatch.nearest <- function(nwps,dmat,method.psorder)
    {
      boot.data.match <- rbind(data.frame(match=1,ps=nwps),dmat[dmat$match==0,])
      bootmatching    <- MatchIt::matchit(match~.,
                                          data     = boot.data.match,
                                          method   = "nearest",
                                          distance = boot.data.match[,"ps"],
                                          m.order  = method.psorder)
      subjid.EC  <- as.numeric(bootmatching$match.matrix)
      return(subjid.EC)
    }

    medps <- function(invet.ps,n.EC)
    {
      n.C   <- length(invet.ps)
      wk.n1 <- floor(n.C/n.EC)
      wk.n2 <- n.C-wk.n1*n.EC
      wk.n3 <- sample(n.EC,wk.n2)
      wk.n4 <- data.frame(X1=1:n.EC,X2=wk.n1)
      wk.n4[wk.n3,"X2"] <- wk.n1+1
      wk.n5 <- wk.n4 %>% mutate(en=cumsum(X2))
      wk.n6 <- data.frame(wk.n5,st=c(1,wk.n5$en[-n.EC]+1))

      wk.ps1 <- sort(invet.ps)
      wk.ps2 <- apply(as.matrix(1:n.EC),1,function(x){median(wk.ps1[wk.n6$st[x]:wk.n6$en[x]])})

      return(wk.ps2)
    }

    ncov     <- length(cov.lab)
    data.cov <- data.frame(data[,c("study","treat",cov.lab)],
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
    data.match <- data.cov.ps[match.pts,c(method.whomatch,"ps")]
    colnames(data.match) <- c("match","ps")

    if(method.whomatch=="conc.contl"){

      n.CC  <- sum(CC.flg)
      ratio <- n.EC/n.CC

      if(method.matching=="optimal"){

        fit.match <- MatchIt::matchit(match~.,
                                      data     = data.match,
                                      method   = method.matching,
                                      distance = data.match[,"ps"],
                                      ratio    = ratio)
        subjid.EC <- as.numeric(fit.match$match.matrix)

      }else if(method.matching=="nearest"){

        fit.match <- MatchIt::matchit(match~.,
                                      data     = data.match,
                                      method   = method.matching,
                                      distance = data.match[,"ps"],
                                      m.order  = method.psorder,
                                      ratio    = ratio)
        subjid.EC <- as.numeric(fit.match$match.matrix)
      }

    }else if((method.whomatch=="conc.treat")|(method.whomatch=="conc.all")){

      if(method.matching=="optimal"){

        fit.match <- MatchIt::matchit(match~.,
                                      data     = data.match,
                                      method   = method.matching,
                                      distance = data.match[,"ps"])
        subjid.EC <- sample(as.numeric(fit.match$match.matrix),n.EC)

      }else if(method.matching=="nearest"){

        ps.M1 <- data.match[data.match[,"match"]==1,"ps"]
        ps.M0 <- data.match[data.match[,"match"]==0,"ps"]
        ps.M  <- sort(abs(apply(as.matrix(ps.M1),1,function(x){x-ps.M0})))

        ncal <- n.EC
        while(1){
          fit.match <- MatchIt::matchit(match~.,
                                        data        = data.match,
                                        method      = method.matching,
                                        distance    = data.match[,"ps"],
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

      }else if((method.matching=="boot.random.optimal")|(method.matching=="boot.random.nearest")|
               (method.matching=="boot.highprob.optimal")|(method.matching=="boot.highprob.nearest")){

        pre.new.ps <- data.match[data.match$match==1,"ps"]
        boot       <- replicate(n.boot,sample(pre.new.ps,n.EC))

        if((method.matching=="boot.random.optimal")|(method.matching=="boot.highprob.optimal")){
          bs.match <- apply(boot,2,function(x){return(bootmatch.optimal(nwps=x,dmat=data.match))})
        }else if((method.matching=="boot.random.nearest")|(method.matching=="boot.highprob.nearest")){
          bs.match <- apply(boot,2,function(x){return(bootmatch.nearest(nwps=x,dmat=data.match,method.psorder=method.psorder))})
        }

        ECp.lab <- as.numeric(rownames(data.match[data.match$match==0,]))
        bc.prob <- apply(as.matrix(ECp.lab),1,function(x){return(mean(bs.match==x))})

        if((method.matching=="boot.random.optimal")|(method.matching=="boot.random.nearest")){
          subjid.EC <- sample(ECp.lab,n.EC,prob=bc.prob)
        }else if((method.matching=="boot.highprob.optimal")|(method.matching=="boot.highprob.nearest")){
          subjid.EC <- ECp.lab[order(bc.prob,decreasing=TRUE)[1:n.EC]]
        }

      }else if((method.matching=="med.optimal")|(method.matching=="med.nearest")|
               (method.matching=="km.optimal")|(method.matching=="km.nearest")|
               (method.matching=="cm.optimal")|(method.matching=="cm.nearest")){

        pre.new.ps <- data.match[data.match$match==1,"ps"]

        if((method.matching=="med.optimal")|(method.matching=="med.nearest")){
          new.ps <- medps(invet.ps=pre.new.ps,n.EC=n.EC)
        }else if((method.matching=="km.optimal")|(method.matching=="km.nearest")){
          new.ps <- as.vector(kmeans(pre.new.ps,n.EC)$centers)
        }else if((method.matching=="cm.optimal")|(method.matching=="cm.nearest")){
          new.ps <- as.vector(e1071::cmeans(pre.new.ps,n.EC)$centers)
        }
        new.data.match <- rbind(data.frame(match=1,ps=new.ps),data.match[data.match$match==0,])

        if((method.matching=="med.optimal")|(method.matching=="km.optimal")|(method.matching=="cm.optimal")){

          fit.match <- MatchIt::matchit(match~.,
                                        data     = new.data.match,
                                        method   = "optimal",
                                        distance = new.data.match[,"ps"])

        }else if((method.matching=="med.nearest")|(method.matching=="km.nearest")|(method.matching=="cm.nearest")){

          fit.match <- MatchIt::matchit(match~.,
                                        data     = new.data.match,
                                        method   = "nearest",
                                        distance = new.data.match[,"ps"],
                                        m.order  = method.psorder)
        }
        subjid.EC <- as.numeric(fit.match$match.matrix)
      }

    }else if(method.whomatch=="treat2contl"){

      if(method.matching=="optimal"){

        fit.match <- MatchIt::matchit(match~.,
                                      data     = data.match,
                                      method   = method.matching,
                                      distance = data.match[,"ps"])
        subjid.C1 <- as.numeric(fit.match$match.matrix)
        subjid.C2 <- subjid.C1[data.cov.ps[subjid.C1,"study"]==0]
        if(length(subjid.C2)>=n.EC){
          subjid.EC <- sample(subjid.C2,n.EC)
        }else{
          subjid.EC <- subjid.C2
        }

      }else if(method.matching=="nearest"){

        ps.M1 <- data.match[data.match[,"match"]==1,"ps"]
        ps.M0 <- data.match[data.match[,"match"]==0,"ps"]
        ps.M  <- sort(abs(apply(as.matrix(ps.M1),1,function(x){x-ps.M0})))

        ncal <- n.EC
        while(1){
          fit.match <- MatchIt::matchit(match~.,
                                        data        = data.match,
                                        method      = method.matching,
                                        distance    = data.match[,"ps"],
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

    return(list(subjid.EC=subjid.EC,data.ps=data.cov.ps[,c("study","treat",cov.lab,"ps")]))
  }
}


#' Summarizing simulation study results
#'
#' Simulation study results are summarized.
#' @usage
#' psborrow.summary(object)
#' @param object List of simulation results.
#' @import stats
#' @rawNamespace import(dplyr,except=c(lag,filter))
#' @export
psborrow.summary <- function(object)
{
  reject  <- object$reject
  theta   <- object$theta
  colflg  <- object$method.borrow
  t.theta <- object$true.theta

  rrate    <- apply(reject[,colflg],2,mean)
  bias     <- apply(theta[theta$measure=="mean",colflg]-t.theta,2,mean)
  empsd    <- apply(theta[theta$measure=="mean",colflg],2,stats::sd)
  modsd    <- apply(theta[theta$measure=="sd",colflg],2,mean)
  empmodsd <- empsd/modsd
  covprob  <- apply((theta[theta$measure=="lcri",colflg]<t.theta)*(theta[theta$measure=="ucri",colflg]>t.theta),2,mean)

  res.out <- data.frame(measure=c("Reject","Bias","EmpSD","ModSD","EmpSD/ModSD","Coverage Prob."),
                        round(rbind(rrate,bias,empsd,modsd,empmodsd,covprob),digits=3))

  colnames(res.out) <- c("measure",colflg)
  rownames(res.out) <- NULL

  ov <- object$ov

  sum.ov <- ov %>%
            group_by(factor,type) %>%
            summarise_at(vars(ov),c(mean,median,sd),na.rm=TRUE)

  res.cov <- data.frame(
               factor = sum.ov$factor,
               type   = sum.ov$type,
               mean   = round(sum.ov$fn1,digits=3),
               median = round(sum.ov$fn2,digits=3),
               sd     = round(sum.ov$fn3,digits=3))

  result <- list(res.out         = res.out,
                 res.cov         = res.cov,
                 n.CT            = object$n.CT,
                 n.CC            = object$n.CC,
                 n.ECp           = object$n.ECp,
                 n.EC            = object$n.EC,
                 drift           = object$drift,
                 method.psest    = object$method.psest,
                 method.pslink   = object$method.pslink,
                 method.whomatch = object$method.whomatch,
                 method.matching = object$method.matching,
                 method.psorder  = object$method.psorder)

  class(result) <- "psborrow.summary"
  result
}

#' @export
print.psborrow.summary <- function(x,...)
{
  ptab1  <- as.matrix(x$res.out[,-1])
  rnames <- c("Reject","Bias","EmpSD","ModSD","EmpSD/ModSD","Coverage Prob.")
  dimnames(ptab1) <- list(rnames,names(x$res.out[,-1]))

  ptab2  <- as.matrix(x$res.cov[,-1])
  rnames <- x$res.cov[,1]
  dimnames(ptab2) <- list(rnames,names(x$res.cov[,-1]))

  ptab3 <- rbind(
    x$n.CT,
    x$n.CC,
    x$n.ECp,
    x$n.EC,
    x$drift)
  rnames <- c("No. of conc. treatment",
              "No. of conc. control",
              "No. of ext. control pool",
              "No. of ext. control",
              "Prior data conflict")
  dimnames(ptab3) <- list(rnames,"")

  if(x$method.matching=="optimal"){
    ptab4 <- rbind(
      x$method.psest,
      x$method.pslink,
      x$method.whomatch,
      x$method.matching)
    rnames <- c("PS estimation",
                "Link function",
                "Who to match",
                "Matching method")
    dimnames(ptab4) <- list(rnames,"")
  }else{
    ptab4 <- rbind(
      x$method.psest,
      x$method.pslink,
      x$method.whomatch,
      x$method.matching,
      x$method.psorder)
    rnames <- c("PS estimation",
                "Link function",
                "Who to match",
                "Matching method",
                "Matching order")
    dimnames(ptab4) <- list(rnames,"")
  }

  cat("\n")
  cat("Treatment effect estimation:\n")
  print.default(ptab1,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Overlapping coefficient:\n")
  print.default(ptab2,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Trial design settings:\n")
  print.default(ptab3,print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Method:\n")
  print.default(ptab4,print.gap=2L,quote=FALSE)
  cat("\n")
}

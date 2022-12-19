
#' Summarizing simulation study results
#'
#' Simulation study results are summarized.
#' @usage
#' commensurate.summary(simres,method.lab,t.theta)
#' @param simres List of simulation results, including reject and theta.
#' @param method.lab Label of information borrowing method.
#' @param t.theta True value of treatment effect.
#' @export

commensurate.summary <- function(simres,method.lab,t.theta)
{
  nmethod <- ncol(simres$reject)

  res.reject <- NULL
  res.bias   <- NULL
  for(i in 1:nmethod){
    res.reject      <- c(res.reject,     round(mean(simres$reject[,i]),digits=3))
    res.mean.bias   <- c(res.mean.bias,  round(mean(simres$theta[,i,1]-t.theta),digits=3))
    res.median.bias <- c(res.median.bias,round(mean(simres$theta[,i,2]-t.theta),digits=3))
    res.emp.sd      <- c(res.emp.sd,     round(sd(simres$theta[,i,1]),digits=3))
    res.mod.sd      <- c(res.mod.sd,     round(mean(simres$theta[,i,3]),digits=3))
  }
  res <- data.frame(X0=c("Reject","MeanBias","MedianBias","EmpSD","ModSD"),
                    rbind(res.reject,res.mean.bias,res.median.bias,res.emp.sd,res.mod.sd))

  colnames(res) <- c("measure",method.lab)
  rownames(res) <- NULL

  return(res)
}

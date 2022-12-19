
#' Summarizing simulation study results
#'
#' Simulation study results are summarized.
#' @usage
#' commensurate.summary(simres,method.lab,t.theta)
#' @param simres List of simulation results, including reject and theta.
#' @param method.lab List of information borrowing method.
#' @param t.theta True value of treatment effect.
#' @export

commensurate.summary <- function(simres,method.lab,t.theta)
{
  nmethod <- ncol(simres$reject)

  res.reject <- NULL
  res.bias   <- NULL
  for(i in 1:nmethod){
    res.reject <- c(res.reject,mean(simres$reject[,i]))
    res.bias   <- c(res.bias,mean(simres$theta[,i,1]-t.theta))
  }
  res <- data.frame(X0=c("reject","bias"),rbind(res.reject,res.bias))
  colnames(res) <- c("measure",method.lab)

  return(res)
}

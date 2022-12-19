% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commensurate.bin.R
\name{commensurate.bin}
\alias{commensurate.bin}
\title{Simulation study of hybrid control design with commensurate power prior:
binary outcome}
\usage{
commensurate.bin(
  n.CT, n.CC, n.EC,
  out.prob.CT, out.prob.CC, driftOR,
  cov.CT, cov.CC, cov.EC, cormat, method,
  chains=2, iter=4000, warmup=floor(iter/2), thin=1,
  alternative="greater", sig.level=0.025, nsim=10)
}
\arguments{
\item{n.CT}{Number of patients in concurrent treatment.}

\item{n.CC}{Number of patients in concurrent control.}

\item{n.EC}{Number of patients in external control.}

\item{out.prob.CT}{True probability of outcome in concurrent treatment.}

\item{out.prob.CC}{True probability of outcome in concurrent control.}

\item{driftOR}{OR between external control and concurrent control for which
the bias should be plotted.}

\item{cov.CT}{List of covariates in concurrent treatment. Distribution and
and its parameters need to be specified for each covariate.}

\item{cov.CC}{List of covariates in concurrent control. Distribution and
and its parameters need to be specified for each covariate.}

\item{cov.EC}{List of covariates in historical control. Distribution and
and its parameters need to be specified for each covariate.}

\item{cormat}{Matrix of correlation coefficients for outcome and covariates,
specified as Gaussian copula parameter.}

\item{method}{List of information borrowing method. \code{"noborrow"} uses
the concurrent data only. \code{"fullborrow"} uses the external control data
without discounting. \code{"cauchy"} uses the commensurate prior to
dynamically borrow the external control data, and the commensurate parameter
is assumed to follow half cauchy distribution. \code{"normal"} uses the
commensurate prior to dynamically borrow the external control data, and the
commensurate parameter is assumed to follow half normal distribution.
\code{"cauchy"} and \code{"normal"} require to specify the scale parameter of
half cauchy and half normal distribution respectively.}

\item{chains}{Number of Markov chains in MCMC sampling. The default value is
\code{chains=2}.}

\item{iter}{Number of iterations for each chain (including warmup) in MCMC
sampling. The default value is \code{iter=4000}.}

\item{warmup}{Number of warmup (aka burnin) iterations per chain in MCMC
sampling. The default value is \code{warmup=floor(iter/2)}.}

\item{thin}{Period for saving samples in MCMC sampling. The default value
is \code{thin=1}.}

\item{alternative}{Alternative hypothesis to be tested ("greater" or "less").
The default value is \code{alternative="greater"}.}

\item{sig.level}{Significance level. The default value is
\code{sig.level=0.025}.}

\item{nsim}{Number of simulated trials. The default value is \code{nsim=10}.}
}
\value{
\item{reject}{\code{TRUE} when significant; otherwise \code{FALSE}.}
\item{theta}{Posterior mean, median, and sd of log odds ratio.}
}
\description{
Bayesian dynamic information borrowing with commensurate power prior is
implemented for hybrid control design, where the concurrent control is
augmented by historical control. No borrowing and full borrowing are also
implemented. The binary outcome is applicable.
}
\examples{
n.CT  <- 100
n.CC  <- 100
n.EC  <- 100

out.prob.CT <- 0.2
out.prob.CC <- 0.2
driftOR     <- 1.0

cov.CT <- list(list(dist="norm", mean=0,sd=1),
               list(dist="binom",prob=0.4))

cov.CC <- list(list(dist="norm", mean=0,sd=1),
               list(dist="binom",prob=0.4))

cov.EC <- list(list(dist="norm", mean=0,sd=1),
            list(dist="binom",prob=0.4))

cormat <- rbind(c(  1,0.1,0.1),
                c(0.1,  1,0.1),
                c(0.1,0.1,  1))

method <- list(list(prior="noborrow"),
               list(prior="fullborrow"),
               list(prior="cauchy",scale=0.5),
               list(prior="normal",scale=0.5))

commensurate.bin(
  n.CT=n.CT, n.CC=n.CC, n.EC=n.EC,
  out.prob.CT=out.prob.CT, out.prob.CC=out.prob.CC, driftOR=driftOR,
  cov.CT=cov.CT, cov.CC=cov.CC, cov.EC=cov.EC, cormat=cormat,
  method=method)
}
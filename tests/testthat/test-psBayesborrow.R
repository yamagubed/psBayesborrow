
n.CT       <- 100
n.CC       <- 50
nevent.C   <- 100
n.ECp      <- 200
nevent.ECp <- 180
n.EC       <- 50
accrual    <- 16

out.mevent.CT <- 6
out.mevent.CC <- 6
driftHR       <- 1

cov.C <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
              list(dist="binom",prob=0.4,lab="cov2"))

cov.cor.C <- rbind(c(  1,0.1),
                   c(0.1,  1))

cov.effect.C <- c(0.9,0.9)

cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
               list(dist="binom",prob=0.4,lab="cov2"))

cov.cor.EC <- rbind(c(  1,0.1),
                    c(0.1,  1))

cov.effect.EC <- c(0.9,0.9)

psmatch.cov <- c("cov1","cov2")

method.whomatch <- "conc.treat"
method.matching <- "optimal"
method.psorder  <- NULL

analysis.cov <- c("cov1")

method.borrow <- list(list(prior="noborrow"),
                      list(prior="normal",scale=0.5))

nsim <- 5

res <- psborrow.t2e(
  n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
  n.ECp=n.ECp, nevent.ECp=nevent.ECp, n.EC=n.EC, accrual=accrual,
  out.mevent.CT=out.mevent.CT, out.mevent.CC=out.mevent.CC, driftHR=driftHR,
  cov.C=cov.C, cov.cor.C=cov.cor.C, cov.effect.C=cov.effect.C,
  cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect.EC=cov.effect.EC,
  psmatch.cov=psmatch.cov, method.whomatch=method.whomatch,
  method.matching=method.matching, method.psorder=method.psorder,
  analysis.cov=analysis.cov, method.borrow=method.borrow,
  chains=1, iter=100, nsim=nsim, seed=100)

rj1.sim <- res$reject[,"noborrow"]
rj2.sim <- res$reject[,"normal0.5"]
th1.sim <- ceiling(100*res$theta[,"noborrow"])/100
th2.sim <- ceiling(100*res$theta[,"normal0.5"])/100
ov.sim  <- ceiling(100*res$ov[,"ov"])/100

rj1.exp <- c(FALSE,FALSE,FALSE,FALSE,FALSE)
rj2.exp <- c( TRUE,FALSE,FALSE,FALSE,FALSE)

th1.exp <- c(
  0.40,0.39,0.27,-0.07,0.87,-0.06,-0.07,0.20,-0.40,0.25,-0.25,-0.22,0.17,-0.68,-0.04,-0.02,-0.01,0.19,-0.40,
  0.30,-0.10,-0.09,0.18,-0.44,0.18)

th2.exp <- c(
  0.27,0.26,0.15,0.01,0.50,-0.06,-0.05,0.27,-0.59,0.34,-0.22,-0.20,0.17,-0.56,0.09,-0.06,-0.07,0.16,-0.39,
  0.23,-0.18,-0.14,0.20,-0.54,0.12)

ov.exp <- c(
  0.90,0.85,0.04,0.87,0.76,-0.01,0.97,0.93,0.00,0.87,0.93,-0.03,0.96,0.91,0.06,0.84,0.86,-0.15,0.91,
  0.91,-0.07,0.93,0.89,0.04,0.97,0.92,0.05,0.92,0.95,-0.08)

test_that("Check psBayesborrow simulation results",{
  expect_identical(rj1.sim,rj1.exp)
  expect_identical(rj2.sim,rj2.exp)
  expect_identical(th1.sim,th1.exp)
  expect_identical(th2.sim,th2.exp)
  expect_identical(ov.sim,  ov.exp)
})

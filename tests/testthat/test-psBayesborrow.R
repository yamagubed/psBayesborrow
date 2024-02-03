
n.CT       <- 100
n.CC       <- 50
nevent.C   <- 100
n.ECp      <- 1000
nevent.ECp <- 800
accrual    <- 16

out.mevent.CT <- 6
out.mevent.CC <- 6
driftHR       <- 1

cov.C <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
             list(dist="binom",prob=0.4,lab="cov2"))

cov.cor.C <- rbind(c(  1,0.1),
                  c(0.1,  1))

cov.EC <- list(list(dist="norm",mean=0,sd=1,lab="cov1"),
            list(dist="binom",prob=0.4,lab="cov2"))

cov.cor.EC <- rbind(c(  1,0.1),
                    c(0.1,  1))

cov.effect <- c(0.1,0.1)

indata <- trial.simulation.t2e(
  n.CT=n.CT, n.CC=n.CC, nevent.C=nevent.C,
  n.ECp=n.ECp, nevent.ECp=nevent.ECp, accrual=accrual,
  out.mevent.CT, out.mevent.CC, driftHR,
  cov.C=cov.C, cov.cor.C=cov.cor.C,
  cov.EC=cov.EC, cov.cor.EC=cov.cor.EC, cov.effect=cov.effect, seed=100)

n.EC <- 50

method.whomatch <- "conc.contl"
method.matching <- "nearest"
method.psorder  <- "large"

mypsmatch <- psmatch(
  study~cov1+cov2, data=indata, n.EC=n.EC,
  method.whomatch=method.whomatch, method.matching=method.matching,
  method.psorder=method.psorder)

subjid.EC.run <- mypsmatch$subjid.EC

subjid.EC.exp <- c(
  657,1108,641,235,960,1123,362,418,967,966,1039,945,1065,615,528,692,475,946,
  1136,242,561,860,808,1119,309,583,662,1147,607,695,1128,1033,488,739,176,164,
  1012,388,670,969,857,982,700,1018,451,448,1121,663,208,938)

test_that("Check psmatch results",{
  expect_identical(subjid.EC.run,subjid.EC.exp)
})

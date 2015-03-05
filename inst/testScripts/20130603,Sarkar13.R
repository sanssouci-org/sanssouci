F <- function(u, v) {
  u*v
}

source("package/R/Sarkar13.R")

n <- 40
gamma <- 0.2
n0 <- 90
K <- 50
k <- 1

beta <- 0.3
arg(n0, n, K, k, beta, gamma, F)

system.time(betaStar <- betaSU(0.1, k, gamma, F, n))

tau <- alphaPrime(1:n, betaStar$root, gamma)

## A quicker mStar?
ms <- sapply(1:n, mStar, gamma)

n <- 1e4;
system.time(ms1 <- sapply(1:n, mStar, gamma))
system.time(ms2 <- mStar2(1:n, gamma))
system.time(ms3 <- mStar3(1:n, gamma))

max(abs(ms2-ms1))
max(abs(ms3-ms1))


n <- 1e6
system.time(ms2 <- mStar2(1:n, gamma))
system.time(ms3 <- mStar3(1:n, gamma))

max(abs(ms3-ms2))

##
Rprof("betaSU.Rout")
betaStar <- betaSU(0.1, k, gamma, F)
Rprof(NULL)
summaryRprof("betaSU.Rout")

m <- 10
pi0 <- 0.6

rho <- 0
p <- 0.5
n <- 20

## nulls and alternatives
m0 <- round(m*pi0)
m1 <- m-m0
H0 <- 1:m0
H1 <- (m0+1):m

## response
y <- rbinom(n, 1, p)

## equi-correlated
if (TRUE) {  ## work fine in fact
  sim <- simulateGaussianNullsFromFactorModel(m, B=n, flavor="equi-correlated", rho=rho)
  str(sim)
  ## image(sim$Sigma)
  eps <- sim$Y
  str(eps)
} else { ## a la mano
  h <- 1
  K <- length(h) ## one factor
  P <- matrix(1, m, K)/sqrt(m) ## K=1 columns
  Z <- matrix(rnorm(m*n), nrow=m, ncol=n)
  W <- matrix(sqrt(h)*rnorm(K), nrow=K, ncol=n)  ## the same 'W' for all 'n' observations !!!
  eps <- sqrt(1-rho)*Z + sqrt(rho) * P %*% W  
}

## means
mu <- matrix(0, nrow=nrow(eps), ncol=ncol(eps)) ## m x n
w1 <- which(y==1)
mu[H1, w1] <- 1
image(mu)

## X
X <- mu+eps

if (FALSE) { ## done directly by package "howmany"
  ## test statistics
  testStat <- function(x, y) {
    ## wilcox.test(x[y==0], x[y==1])$statistic
    wilcox.test(x[y==0], x[y==1])$p.value
  }
  stats <- apply(X, 1, testStat, y)
  plot(stats)
  
  ## permutation
  yB <- sample(y)
  statsB <- apply(X, 1, testStat, yB)

  ## permutations
  B <- 500
  statsB <- replicate(B, {
    yB <- sample(y)
    statsB <- apply(X, 1, testStat, yB)
  })
  dim(statsB)
}

library(howmany.pn)
Rprof("howmany.Rout")
set.seed(123)
res <- howmany_dependent(t(X), y)
Rprof(NULL)
summaryRprof("howmany.Rout")

plot(res)
o <- res$order
## true discoveries
S <- cumsum(o>m0)
Sbar <- lowerbound(res)
plot(S, Sbar, t='s')
abline(a=0, b=1)

## faster ??
source("package/inst/testScripts/mywilcox.test.R")

Rprof("howmany,mine.Rout")
set.seed(123)
resM <- howmany_dependent(t(X), y, test=mywilcox.test)
Rprof(NULL)
summaryRprof("howmany,mine.Rout")$sampling.time

Rprof("howmany,mine2.Rout")
set.seed(123)
resM2 <- howmany_dependent(t(X), y, test=function(x, y, ...) mywilcox.test(x, y, p.value=FALSE))
Rprof(NULL)
summaryRprof("howmany,mine2.Rout")$sampling.time

Rprof("howmany,mine3.Rout")
set.seed(123)
resM3 <- howmany_dependent(t(X), y, test=function(x, y, ...) mywilcox.test(x, y, p.value=FALSE, bypassTiesMethod=TRUE))
Rprof(NULL)
summaryRprof("howmany,mine3.Rout")$sampling.time

summaryRprof("howmany,mine.Rout")$sampling.time/summaryRprof("howmany.Rout")$sampling.time
summaryRprof("howmany,mine2.Rout")$sampling.time/summaryRprof("howmany.Rout")$sampling.time
summaryRprof("howmany,mine3.Rout")$sampling.time/summaryRprof("howmany.Rout")$sampling.time
## ~20% of the original time.

## same result ?
resMM <- resM2

om <- resMM$order
Sm <- cumsum(om>m0)
Sbarm <- lowerbound(resMM)


identical(om, o)  ## not always true because several (quite null) hypotheses have a test statistic of 0 exactly...
identical(Sm, S)
identical(Sbar, Sbarm)

## next step would require updating the package "howmany" to calculate the permuted test statistics more efficiently (sorting X just once and accessing its elements)

library(sansSouci)
source(system.file("testScripts/Mein2006/R/wilcox.stat.R", package="sansSouci"))

foo <- function(m, n, B, alpha) {
  ## independence
  sim <- simulateGaussianNullsFromFactorModel(m, B=n, flavor="equi-correlated", rho=0, cov=FALSE)
  eps <- sim$Y

  ## means under the total null: 0
  mu <- matrix(0, nrow=nrow(eps), ncol=ncol(eps)) ## m x n
  X <- mu+eps

  y <- rbinom(n, 1, 1/2)
  
  w <- wilcox.stat(X, y, B=B)
  scoreMat <- w$stat0Mat
  o <- order(w$stat, decreasing=TRUE)
  stat <- w$stat[o]

  ## joint FWER control (through gammatification)
  resJ <- getJointFWERThresholds(scoreMat, tau="kFWER", alpha=alpha, maxSteps=1000)
  ## str(resJ)
  thr <- resJ$thr
  gg <- sum(stat>thr)
  ge <- sum(stat>=thr)
  prob <- resJ$prob

  BB <- sapply(stat, function(x) sum(x<=thr))  ## Eqn (7) in Meinshausen (2006) (*not* Vbar)
  R <- 1:m
  Sbar <- pmax(0, cummax(R-BB[R]))

  BB1 <- sapply(stat, function(x) sum(x<thr))  ## Eqn (7) in Meinshausen (2006) (*not* Vbar)
  Sbar1 <- pmax(0, cummax(R-BB1[R]))

  c(gg, ge, prob, Sbar=sum(Sbar>0), Sbar1=sum(Sbar1>0))
}

m <- 1e3
n <- 1e2
B <- 1e4
alpha <- 0.2

nbCores <- 10
nbSimu <- 2*nbCores
system.time(res <- parallel::mclapply(1:nbSimu, FUN=function(ii) {
  print(ii);
  foo(m, n, B, alpha)
}, mc.cores=nbCores))
res <- t(simplify2array(res))
head(res)

colSums(res[, -3]>0)/nrow(res)
summary(res[, 3])
stem(res[, 3]) ## seems OK for B large enough

head(res[, c(1, 4)], 20)
head(res[, c(2, 5)], 20)

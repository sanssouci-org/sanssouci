simulateMein2006 <- structure(function(
    m,
    ### Number of hypotheses
    rho,
    ### Level of equi-correlation between pairs of variables
    n,
    ### Number of observations
    pi0,
    ### Proportion of true null hypotheses
    SNR=1,
    ### Signal to noise ratio. Either a numeric value (a measure of distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
    p=0.5,
    ### Probability of success of the outcome variable
    w=NULL
### An optional vector of length \code{n}, the underlying factor
### driving equi-correlation
    ) {
  m0 <- round(m*pi0)
  m1 <- m-m0
  H <- rep(c(0, 1), times=c(m0, m1))
  H0 <- which(H==0)
  H1 <- which(H==1)
  if (length(SNR)>1) {
    stopifnot(length(SNR)==m1)
  }
  
  ## equi-correlated noise
  eps <- simulateGaussianEquiCorrelatedNulls(m, n=n, rho=rho, w=w)
  w <- attr(eps, "w")
  
  ## binomial response
  y <- rbinom(n, 1, p)

  ## means
  mu <- matrix(0, nrow=nrow(eps), ncol=ncol(eps)) ## m x n
  if (m0<m) {
    w1 <- which(y==1)
    mu[H1, w1] <- SNR*sqrt(2*log(n)/n)
  }
  X <- mu+eps
  list(
      X=X,
### An \eqn{m x n} covariate matrix
      y=y,
### A vector of \eqn{n} phenotypes in \eqn{{0,1}}
      H=H,
### A vector of length \eqn{m}, the status of each
### hypothesis:\describe{
### \item{0}{true null hypothesis}
### \item{1}{true alternative hypothesis}
### }
      w=w)
### A vector of length \code{n}, the underlying factor
### driving equi-correlation
}, ex=function() {
  m <- 123
  rho <- 0.2
  n <- 100
  pi0 <- 0.5

  sim <- simulateMein2006(m, rho, n, pi0, SNR=1)
  X <- sim$X
  y <- sim$y
  
  w <- wilcoxStat(X, y, B=B)
  scoreMat <- w$stat0Mat
  stat <- w$stat

  ## show test statistics
  pch <- 20
  colStat <- 1+sim$H
  plot(stat, col=colStat, main="Test statistics", pch=pch)
  legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
})

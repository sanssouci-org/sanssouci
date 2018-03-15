##' Simulate test statistics as in Meinshausen (2006)
##'
##' @param m Number of hypotheses
##' @param rho Level of equi-correlation between pairs of variables
##' @param n Number of observations
##' @param B Number of resamplings to estimate the test statistics
##' @param pi0 Proportion of true null hypotheses
##' @param SNR Signal to noise ratio. Either a numeric value (a measure of
##' distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
##' @param p Probability of success of the outcome variable
##' @param w An optional vector of length \code{n}, the underlying factor
##' driving equi-correlation
##' @return A list with elements \describe{
##' \item{x}{A vector of length \eqn{m} test statistics}
##' \item{X0}{An \eqn{m x B} matrix of test statistics under the null
##' hypothesis}
##' \item{H}{A vector of length \eqn{m}, the status of each
##' hypothesis: 0 for true null hypothesis, and 1 for true alternative
##' hypothesis} }
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @importFrom stats rbinom
##' @examples
##'
##' m <- 123
##' rho <- 0.2
##' n <- 100
##' pi0 <- 0.5
##' B <- 1e3
##'
##' sim <- sansSouci:::simulateMein2006(m, rho, n, B, pi0, SNR=1)
##' scoreMat <- sim$X0
##' stat <- sim$x
##'
##' ## show test statistics
##' pch <- 20
##' colStat <- 1+sim$H
##' plot(stat, col=colStat, main="Test statistics", pch=pch)
##' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
##'
simulateMein2006 <- structure(function(m, rho, n, B, pi0, SNR=1, p=0.5, w=NULL) {
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
    ## *not* faster:
    ## mu <- SNR*sqrt(2*log(n)/n)
    ## if (m0<m) {
    ##     w1 <- which(y==1)
    ##     eps[H1, w1] <- eps[H1, w1] + mu
    ## }
    
    w <- wilcoxStat(X, y, B=B)
    X0 <- w$stat0Mat
    x <- w$stat
    list(x=x,
         X0=X0,
         H=H)
}, ex=function() {
    m <- 123
    rho <- 0.2
    n <- 100
    pi0 <- 0.5
    B <- 1000
    
    sim <- simulateMein2006(m, rho, n, B, pi0, SNR=1)
    scoreMat <- sim$X0
    stat <- sim$x
    
    ## show test statistics
    pch <- 20
    colStat <- 1+sim$H
    plot(stat, col=colStat, main="Test statistics", pch=pch)
    legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
})

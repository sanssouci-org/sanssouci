##' Simulate Gaussian equi-correlated test statistics
##'
##' @param m Number of hypotheses
##' @param rho Level of equi-correlation between pairs of variables
##' @param B Number of simulations
##' @param pi0 Proportion of true null hypotheses
##' @param SNR Signal to noise ratio. Either a numeric value (a measure of
##' distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
##' @return A list with elements \describe{
##' \item{x}{A vector of length \eqn{m} test statistics}
##' \item{X0}{An \eqn{m x B} matrix of test statistics under the null
##' hypothesis}
##' \item{H}{A vector of length \eqn{m}, the status of each
##' hypothesis: 0 for true null hypothesis, and 1 for true alternative
##' hypothesis} }
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @examples
##'
##' m <- 123
##' rho <- 0.2
##' B <- 100
##' pi0 <- 0.5
##'
##' sim <- simulateEqui(m, rho, B, pi0, SNR=1)
##' stat <- sim$x
##'
##' ## show test statistics
##' pch <- 20
##' colStat <- 1+sim$H
##' plot(stat, col=colStat, main="Test statistics", pch=pch)
##' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
##'
simulateEqui <- function(m, rho, B, pi0, SNR=1){
    m0 <- round(m*pi0)
    m1 <- m-m0
    H <- rep(c(0, 1), times=c(m0, m1))
    H <- sample(H)
    H1 <- which(H==1)
    if (length(SNR)>1) {
        stopifnot(length(SNR)==m1)
    }
    ## signals
    mu <- rep(0, m)
    mu[H1] <- SNR

    ## equi-correlated noise
    sim <- simulateGaussianEquiCorrelatedNulls(m, n=1+B, rho=rho)
    x <- mu + sim[, 1]
    Xb <- sim[, -1, drop=FALSE]
    list(x=x,
         X0=Xb,
         H=H)
}

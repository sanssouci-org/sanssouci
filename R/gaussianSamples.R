#' Simulate equi-correlated data
#' 
#' Simulate equi-correlated data coming from one or two populations
#'
#' @param m Number of hypotheses
#' @param rho Level of equi-correlation between pairs of variables
#' @param n Number of observations, i.e. sample size
#' @param pi0 Proportion of true null hypotheses
#' @param SNR Signal to noise ratio. Either a numeric value (a measure of
#'   distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
#' @param p A numeric value in (0, 1], the frequency of the population under H1. If \eqn{p=1} (the default), a data set from a single population is generated.
#' @param w An optional vector of length \code{n}, the underlying factor driving
#'   equi-correlation
#' @return A list with two elements: \describe{
#' \item{X}{An \eqn{m x n} matrix of \eqn{m}-dimensional Gaussian observations, where \eqn{m1} rows are under H1, and \eqn{m0} rows are under H0, with \eqn{mO/m = pi0}}
#' \item{H}{A vector of length \eqn{m}, the status of each
#' hypothesis: 0 for true null hypothesis, and 1 for true alternative
#' hypothesis}}.
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @details. If \code{p = 1}, each of m1 variables under H1 has mean \eqn{SNR/\sqrt(n)}. If \code{0 < p < 1} n0 samples are drawn from a N(0,1) distribution while n1 are drawn from a N(mu, 1) distribution, where \eqn{mu = SNR*sqrt(1/(n0) + 1/n1)}. The argument \code{p} is the probability of a sample to belong to the non-zero-mean population
#' @export
#' @importFrom stats rbinom
#' @examples
#'
#' m <- 123
#' rho <- 0.2
#' n <- 100
#' pi0 <- 0.5
#' B <- 1e3
#'
#' ## two-sample data
#' sim <- gaussianSamples(m, rho, n, pi0, SNR=2, p=0.5)
#' tests <- testByRandomization(sim$X, B, flavor = "perm")
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col=colStat, main="Test statistics", pch=pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#'
#' ## sign-flipping
#' sim <- gaussianSamples(m, rho, n, pi0, SNR=2)
#' tests <- testByRandomization(sim$X, B, flavor = "flip")
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col=colStat, main="Test statistics", pch=pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#'
gaussianSamples <- function(m, rho, n, pi0, SNR = 1, p = 1, w = NULL) {
    m0 <- round(m*pi0)
    m1 <- m - m0
    H <- rep(c(0, 1), times = c(m0, m1))
    H1 <- which(H == 1)

    ## sanity checks
    if (length(SNR) > 1) {
        stopifnot(length(SNR)==m1)
    }

    ## 1. equi-correlated noise
    eps <- simulateGaussianEquiCorrelatedNulls(m, n = n, rho = rho, w = w)
    w <- attr(eps, "w")

    ## 2. signal
    mu <- matrix(0, nrow = nrow(eps), ncol = ncol(eps)) ## m x n
    if (m0 < m) {
        if (p == 1) {  
            # one-sample data (defined here as a corner case of the two-sample problem)
            signal <- SNR/sqrt(n)
            mu[H1, ] <- signal
            colnames(mu) <- rep(1, n)
        } else {
            y <- rbinom(n, 1, p)     ## Bernoulli response
            w1 <- which(y == 1)
            n1 <- length(w1)
            signal <- SNR*sqrt(1/(n-n1) + 1/n1)  ## ! scaling is test-specific (here T)!
            mu[H1, w1] <- signal
            colnames(mu) <- y
        }
    }

    # data: signal + noise
    X <- mu + eps
    list(X = X, H = H)
}



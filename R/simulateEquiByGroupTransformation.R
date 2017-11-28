#' Simulate test statistics by group transformation
#' 
#' @param m Number of hypotheses
#' @param rho Level of equi-correlation between pairs of variables
#' @param n Number of observations
#' @param B Number of resamplings to estimate the test statistics
#' @param pi0 Proportion of true null hypotheses
#' @param flavor A character value, the type of group transformation to be
#'   performed. Should either be "perm" for two-sample permutation or "flip" for
#'   sign flipping
#' @param SNR Signal to noise ratio. Either a numeric value (a measure of 
#'   distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
#' @param p Probability of success of the outcome variable for flavor 
#'   "perm"
#' @param w An optional vector of length \code{n}, the underlying factor driving
#'   equi-correlation
#' @details
#' 
#' For flavor "two-sample", we test the null hypothesis: "both groups have the 
#' same mean" against the one-sided alternative that the mean is larger in the 
#' second group. We use a Student test statistic, although other statistics such
#' as the Mann-Whitney statistic could be used as well. Permuted test statistics
#' are calculated by B permutations of the group labels. Corresponding observed 
#' and permuted p-values are calculated as the proportion of permutations 
#' (including the identity) for which the permuted test statistic is larger than
#' the observed test statistic.
#' 
#' For flavor "flip", we test the null hypothesis: "the mean is 0" against the 
#' two-sided alternative that the mean is larger than 0. We use the (rescaled) 
#' empirical mean of the observations as a test statistic. Sign-flipped test 
#' statistics are calculated by flipping the sign of each observation with 
#' probability 1/2.
#' 
#' @return A list with elements \describe{ \item{x}{A vector of length \eqn{m} 
#'   test statistics} \item{X0}{An \eqn{m x B} matrix of test statistics under 
#'   the null hypothesis} \item{H}{A vector of length \eqn{m}, the status of 
#'   each hypothesis: 0 for true null hypothesis, and 1 for true alternative 
#'   hypothesis}} The test statistics are \eqn{\sim N(0,1)}, and \eqn{\sim 
#'   N(\mu,1)}, with \eqn{\mu>0}
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
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
#' ## two-sample permutation
#' sim <- simulateEquiByGroupTransformation(m, rho, n, B, pi0, SNR=2, flavor="perm")
#' scoreMat <- sim$X0
#' stat <- sim$x
#' 
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(stat, col=colStat, main="Test statistics", pch=pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#' 
#' ## sign-flipping
#' sim <- simulateEquiByGroupTransformation(m, rho, n, B, pi0, SNR=2, flavor="flip")
#' scoreMat <- sim$X0
#' stat <- sim$x
#' 
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(stat, col=colStat, main="Test statistics", pch=pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#' 
simulateEquiByGroupTransformation <- function(m, rho, n, B, pi0, 
                                           flavor=c("perm", "flip"),
                                           p.value=FALSE,
                                           SNR=1, p=0.5, w=NULL) {
    m0 <- round(m*pi0)
    m1 <- m - m0
    H <- rep(c(0, 1), times=c(m0, m1))
    H1 <- which(H == 1)
    
    ## sanity checks
    if (length(SNR) > 1) {
        stopifnot(length(SNR)==m1)
    }
    flavor <- match.arg(flavor)
    if (flavor == "perm") {
        stopifnot(0 < p && p < 1)        
    }
    ## 1.equi-correlated noise
    eps <- simulateGaussianEquiCorrelatedNulls(m, n = n, rho = rho, w = w)
    w <- attr(eps, "w")
    
    ## 2. signal
    mu <- matrix(0, nrow = nrow(eps), ncol = ncol(eps)) ## m x n
    if (flavor == "perm") {
        if (m0 < m) {
            y <- rbinom(n, 1, p)     ## binomial response
            w1 <- which(y == 1)
            mu[H1, w1] <- SNR*sqrt(2*log(n)/n)
        }
        X <- mu + eps   # data: signal + noise
        tests <- testByTwoSamplePermutation(X, y, B, p.value = p.value, seed=NULL)
    } else if (flavor == "flip") {
        if (m0 < m) {
            mu[H1, ] <- SNR*sqrt(2*log(n)/n)
        }
        X <- mu + eps   # data: signal + noise
        tests <- testBySignFlipping(X, B, p.value = p.value, seed = NULL)
    }
    if (p.value) {      # map to N(0,1), cf issue #2
        res <- list(x = qnorm(1 - tests$p),    
                    X0 = qnorm(1 - tests$p0),
                    H = H)
    } else {
        res <- list(x = tests$T,
                    X0 = tests$T0,
                    H = H)
    }
    res
}

#' Simulate Gaussian test statistics
#'
#' @param m Number of hypotheses
#' @param B Number of simulations
#' @param pi0 Proportion of true null hypotheses
#' @param SNR Signal to noise ratio. Either a numeric value (a measure of
#' distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
#' @param dep A character value, the type of dependency between test statistics. Can be one of "equi" for equi-correlation, or "Toeplitz". Defaults to "equi".
#' @param param A numeric value defaulting to \code{0}. If \code{dep=="equi"}, \code{param} is the level of equi-correlation between pairs of variables. If \code{dep=="Toeplitz"}, the first row of the Toeplitz matrix will be \code{(1:m)^(param)}. 
#' @return A list with elements \describe{
#' \item{x}{A vector of length \eqn{m} test statistics}
#' \item{X0}{An \eqn{m x B} matrix of test statistics under the null
#' hypothesis}
#' \item{H}{A vector of length \eqn{m}, the status of each
#' hypothesis: 0 for true null hypothesis, and 1 for true alternative
#' hypothesis} }
#' @details B
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @importFrom Matrix Matrix
#' @examples
#'
#' m <- 123
#' B <- 100
#'
#' # independent statistics under the full null
#' sim <- gaussianTestStatistics(m, B)
#' 
#' # equi-correlated statistics under the full null
#' sim <- gaussianTestStatistics(m, B, dep = "equi", param = 0.2)
#' 
#' # equi-correlated statistics with some signal
#' sim <- gaussianTestStatistics(m, B, pi0 = 0.8, SNR = 1, dep = "equi", param = 0.2)
#'
#' ## show test statistics
#' stat <- sim$x
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(stat, col=colStat, main="Test statistics", pch=pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#'
#' # Toeplitz statistics with some signal
#' sim <- gaussianTestStatistics(m, B, pi0 = 0.8, SNR = 1, dep = "Toeplitz", param = -0.5)
#' 
gaussianTestStatistics <- function(m, B, pi0 = 1, SNR = 0, dep = c("equi", "Toeplitz"), param = 0){
    # sanity checks
    dep <- match.arg(dep)
    if (dep == "equi") {
        if (param < 0 || param > 1) {
            stop("'param' should be in [0,1] if 'dep' is 'equi'")
        }
    } else if (dep == "Toeplitz") {
        if (param >= 0) {
            stop("'param' should be negative if 'dep' is 'Toeplitz'")
        }
    }
    
    # multiple testing setting
    m0 <- round(m*pi0)
    m1 <- m-m0
    H <- rep(c(0, 1), times = c(m0, m1))
    H <- sample(H)
    H1 <- which(H == 1)

    # signal
    mu <- rep(0, m)
    if (m1 > 0) {
        if (length(SNR) > 1) {
            stopifnot(length(SNR) == m1)
        }
        mu[H1] <- SNR
    }
    
    # noise
    if (dep == "equi") {
        sim <- simulateGaussianEquiCorrelatedNulls(m, n = 1 + B, rho = param)
    } else if (dep == "Toeplitz") {
        tcoefs <- stats::toeplitz((1:m)^(param))
        Sigma <- Matrix(tcoefs, sparse = TRUE)
        sim <- simulateGaussianNullsFromSigma(Sigma, n = 1 + B)
    }
    
    x <- mu + sim[, 1]
    Xb <- sim[, -1, drop = FALSE]
    list(x = x,
         X0 = Xb,
         H = H)
}

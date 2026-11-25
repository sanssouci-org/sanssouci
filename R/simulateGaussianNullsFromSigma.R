#' sim_null_gaussian_sigma
#'
#' Simulate m observations of a p-dimensional Gaussian vector with given
#' covariance matrix
#'
#' @param SIGMA Population covariance matrix (of size \code{m x m})
#' @param n Number of replications of the simulation
#' @return A \code{m x n} \code{Matrix} simulated test statistics, where
#' \code{m} is the dimension of \code{SIGMA}
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom Matrix Matrix isSymmetric
#' @importFrom stats rnorm
#' @keywords internal
#' @export
#' @examples
#'
#' library("Matrix")
#' m <- 100
#' n <- 1000
#'
#' ##' Toeplitz, short range
#' tcoefs <- toeplitz((1:m)^(-2))
#' SIGMA <- Matrix(tcoefs, sparse = TRUE)
#' Y <- sanssouci:::sim_null_gaussian_sigma(SIGMA, n)
#' str(Y)
#'
#' SIGMAHat <- Y %*% t(Y)/n
#' image(SIGMAHat)
#' image(SIGMA)
#'
#' svd(SIGMAHat)$d
#' (ss <- svd(SIGMA)$d)
#' max(ss)/min(ss)
#'
#'
#' ##' Toeplitz, long range
#' tcoefs <- toeplitz((1:m)^(-.1))
#' SIGMA <- Matrix(tcoefs, sparse = TRUE)
#' Y <- sanssouci:::sim_null_gaussian_sigma(SIGMA, n)
#' str(Y)
#'
#' SIGMAHat <- Y %*% t(Y)/n
#' image(SIGMAHat)
#' image(SIGMA)
#'
#' svd(SIGMAHat)$d
#' (sl <- svd(SIGMA)$d)
#' max(sl)/min(sl)
#'
#' plot(ss/sum(ss), col=1, t='b', ylim=c(0,1))
#' lines(sl/sum(sl), col=2, t='b')
#'
sim_null_gaussian_sigma <- function(SIGMA, n = 1) {
  stopifnot(isSymmetric(SIGMA))
  m <- ncol(SIGMA)
  t(chol(SIGMA)) %*% matrix(rnorm(m * n), nrow = m, ncol = n)
}

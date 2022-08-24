# simulateGaussianNullsFromSigma
#
# Simulate m observations of a p-dimensional Gaussian vector with given
# covariance matrix
#
#
# @param Sigma Population covariance matrix (of size \code{m x m})
# @param n Number of replications of the simulation
# @return A \code{m x n} \code{Matrix} simulated test statistics, where
# \code{m} is the dimension of \code{Sigma}
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom Matrix Matrix isSymmetric
#' @importFrom stats rnorm
# @examples
#
# library("Matrix")
# m <- 100
# n <- 1000
#
# ## Toeplitz, short range
# tcoefs <- toeplitz((1:m)^(-2))
# Sigma <- Matrix(tcoefs, sparse = TRUE)
# Y <- sanssouci:::simulateGaussianNullsFromSigma(Sigma, n)
# str(Y)
#
# SigmaHat <- Y %*% t(Y)/n
# image(SigmaHat)
# image(Sigma)
#
# svd(SigmaHat)$d
# (ss <- svd(Sigma)$d)
# max(ss)/min(ss)
#
#
# ## Toeplitz, long range
# tcoefs <- toeplitz((1:m)^(-.1))
# Sigma <- Matrix(tcoefs, sparse = TRUE)
# Y <- sanssouci:::simulateGaussianNullsFromSigma(Sigma, n)
# str(Y)
#
# SigmaHat <- Y %*% t(Y)/n
# image(SigmaHat)
# image(Sigma)
#
# svd(SigmaHat)$d
# (sl <- svd(Sigma)$d)
# max(sl)/min(sl)
#
# plot(ss/sum(ss), col=1, t='b', ylim=c(0,1))
# lines(sl/sum(sl), col=2, t='b')
#
simulateGaussianNullsFromSigma <- function(Sigma, n = 1){
  stopifnot(isSymmetric(Sigma))
  m <- ncol(Sigma)
  ## cs <- chol(Sigma)
  ## replicate(n, cs %*% matrix(rnorm(m), nrow=m, ncol=1))
  t(chol(Sigma)) %*% matrix(rnorm(m*n), nrow = m, ncol = n)
}

###########################################################################
# HISTORY:
# 2013-04-26
# o Created.
###########################################################################


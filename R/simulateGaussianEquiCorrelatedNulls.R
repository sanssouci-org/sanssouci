##' simulateGaussianEquiCorrelatedNulls
##'
##' Simulate equi-correlated null hypotheses
##'
##'
##' @param m Number of tests
##' @param n Number of replications
##' @param rho \code{1-rho} is the standard deviation of the noise
##' @param w An optional vector of length \code{n}, the underlying factor
##' driving equi-correlation
##' @return A \code{m x n} matrix of simulated observations. The underlying
##' factor driving equi-correlation is given as \code{attr(Y, "w")}.
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @examples
##'
##' m <- 1331
##' n <- 111
##' ## equi-correlated
##' rho <- 0.2
##'
##' Y <- sansSouci:::simulateGaussianEquiCorrelatedNulls(m, n, rho)
##' ## check equi-correlation:
##' covmat <- cov(t(Y))
##' diag(covmat) <- NA
##' dim(covmat) <- NULL
##' summary(covmat)
##'
simulateGaussianEquiCorrelatedNulls <- structure(function(
### Simulate equi-correlated null hypotheses
    m,
### Number of tests
    n=1e3,
### Number of replications
    rho=0,
### \code{1-rho} is the standard deviation of the noise
    w=NULL
### An optional vector of length \code{n}, the underlying factor
### driving equi-correlation
    ) {

  Z <- matrix(rnorm(m*n), ncol=n)
  if (rho==0) {
    Y <- Z
  } else {
    if (is.null(w)) {
      w <- rnorm(n)
    } else {
      stopifnot(length(w)==n)
    }
    W <- matrix(w, ncol=n, nrow=m, byrow=TRUE)
    Y <- sqrt(1-rho)*Z + sqrt(rho)*W
  }
  attr(Y, "w") <- w
  Y
### A \code{m x n} matrix of simulated observations. The underlying
### factor driving equi-correlation is given as \code{attr(Y, "w")}.
}, ex=function(){
  m <- 1331
  n <- 111
  ## equi-correlated
  rho <- 0.2

  Y <- simulateGaussianEquiCorrelatedNulls(m, n, rho)
  ## check equi-correlation:
  covmat <- cov(t(Y))
  diag(covmat) <- NA
  dim(covmat) <- NULL
  summary(covmat)
})

############################################################################
## HISTORY:
## 2014-04-21
## o Created from 'simulateFactorModelNullsFromSingularValuesAndLoadings'.
############################################################################


##' simulateFactorModelNullsFromSingularValuesAndLoadings
##'
##' Simulate null hypotheses according to a factor model
##'
##'
##' @param m Number of tests
##' @param h A vector of \code{k} singular values associated to each factor
##' @param P A \code{m x k} matrix of factor loadings
##' @param rho \code{1-rho} is the standard deviation of the noise
##' @return \item{Y}{a vector of \code{m} simulated observations} \item{W}{a
##' vector of \code{k} factors}
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @importFrom Matrix Matrix
##' @export
##' @examples
##'
##' m <- 10
##'
##' ## independent
##' sim0 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m)
##' str(sim0)
##' sum(is.na(simulateFactorModelNullsFromSingularValuesAndLoadings(m)))
##' S0 <- getFactorModelCovarianceMatrix(m)
##' image(S0)
##'
##' ## equi-correlated
##' rho <- 0.2
##' h <- 1
##' P <- matrix(1, m, length(h))
##' sim1 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
##' str(sim1)
##' S1 <- getFactorModelCovarianceMatrix(m, h, P, rho=rho)
##' image(S1)
##'
##' ## 3-factor model
##' m <- 4*floor(m/4) ## make sure m/4 is an integer
##' rho <- 0.5
##' h <- c(0.5, 0.3, 0.2)*m
##' gamma1 <- rep(1,m)
##' gamma2 <- rep(c(-1, 1), each=m/2)
##' gamma3 <- rep(c(-1, 1), times=m/2)
##' P <- cbind(gamma1, gamma2, gamma3)/sqrt(m)
##'
##' sim3 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
##' str(sim3)
##' S3 <- getFactorModelCovarianceMatrix(m, h, P, rho=rho)
##' image(S3)
##' sim3
##'
simulateFactorModelNullsFromSingularValuesAndLoadings <- structure(function(
### Simulate null hypotheses according to a factor model
    m,
### Number of tests
    h=numeric(0),
### A vector of \code{k} singular values associated to each factor
    P=Matrix(nrow=m, ncol=length(h)),
### A \code{m x k} matrix of factor loadings
    rho=0
### \code{1-rho} is the standard deviation of the noise
    ) {
  k <- length(h)
  ## sanity checks
  stopifnot(nrow(P)==m)
  stopifnot(ncol(P)==k)
  if (k>1) {
    ## is P'P orthonormal ?
    mm <- max(abs(diag(k) - t(P) %*% P))
    if (mm > 1e-10) {
      stop("t(P) %*% P should be orthonormal")
    }
    rm(mm)
  }

  Z <- matrix(rnorm(m), ncol=1)
  W <- rnorm(k)
  if (rho==0) {
    Y <- Z
  } else {
    Y <- sqrt(1-rho)*Z + sqrt(rho) * P %*% Matrix(sqrt(h)*W, ncol=1)
#    if (sum(is.na(Y))) browser()
  }
  Y <- as.vector(Y)
  list(Y=Y, W=W)
### \item{Y}{a vector of \code{m} simulated observations}
### \item{W}{a vector of \code{k} factors}
}, ex=function(){
  m <- 10

  ## independent
  sim0 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m)
  str(sim0)
  sum(is.na(simulateFactorModelNullsFromSingularValuesAndLoadings(m)))
  S0 <- getFactorModelCovarianceMatrix(m)
  image(S0)

  ## equi-correlated
  rho <- 0.2
  h <- 1
  P <- matrix(1, m, length(h))
  sim1 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
  str(sim1)
  S1 <- getFactorModelCovarianceMatrix(m, h, P, rho=rho)
  image(S1)

  ## 3-factor model
  m <- 4*floor(m/4) ## make sure m/4 is an integer
  rho <- 0.5
  h <- c(0.5, 0.3, 0.2)*m
  gamma1 <- rep(1,m)
  gamma2 <- rep(c(-1, 1), each=m/2)
  gamma3 <- rep(c(-1, 1), times=m/2)
  P <- cbind(gamma1, gamma2, gamma3)/sqrt(m)

  sim3 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
  str(sim3)
  S3 <- getFactorModelCovarianceMatrix(m, h, P, rho=rho)
  image(S3)
  sim3
})


############################################################################
## HISTORY:
## 2013-03-20
## o Created from Etienne's code.
############################################################################


# getFactorModelCovarianceMatrix
#
# calculate the covariance matrix of a factor model
#
#
# @param m Number of tests
# @param h A vector of \code{k} singular values associated to each factor
# @param P A \code{m x k} matrix of factor loadings
# @param rho \code{1-rho} is the standard deviation of the noise
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom Matrix Matrix
getFactorModelCovarianceMatrix <- function(m, h = numeric(0), 
                                           P = Matrix(nrow=m, ncol=length(h)), rho = 0) {
  k <- length(h)
  ## sanity checks
  stopifnot(nrow(P)==m)
  stopifnot(ncol(P)==k)
  if (k>1) {
    ## is P'P orthonormal ?
    mm <- max(abs(diag(k) - t(P) %*% P))
    if (mm > 1e10) {
      stop("t(P) %*% P should be orthonormal")
    }
    rm(mm)
    mat <-   (1-rho)*diag(1, m)+ rho* P %*% diag(h) %*% t(P)
  } else { ## k==0
    mat <- diag(1, m)
  }
  mat
}

############################################################################
# HISTORY:
# 2014-03-03
# o BUG FIX in the independent case ('P' was not recognized as a matrix).
# 2013-03-20
# o Created from Etienne's code.
############################################################################


#' sim_null_gaussian_equi
#'
#' Simulate equi-correlated null hypotheses
#'
#'
#' @param m Number of tests
#' @param n Number of replications
#' @param rho \code{1-rho} is the standard deviation of the noise
#' @param w An optional vector of length \code{n}, the underlying factor
#' driving equi-correlation
#' @return A \code{m x n} matrix of simulated observations. The underlying
#' factor driving equi-correlation is given as \code{attr(Y, "w")}.
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @keywords internal
#' @examples
#'
#' m <- 1331
#' n <- 111
#' ## equi-correlated
#' rho <- 0.2
#'
#' Y <- sanssouci:::sim_null_gaussian_equi(m, n, rho)
#' ## check equi-correlation:
#' covmat <- cov(t(Y))
#' diag(covmat) <- NA
#' dim(covmat) <- NULL
#' summary(covmat)
#'
sim_null_gaussian_equi <- function(m, n = 1e3, rho = 0, w = NULL) {
  Z <- matrix(rnorm(m * n), ncol = n)
  if (rho == 0) {
    Y <- Z
  } else {
    if (is.null(w)) {
      w <- rnorm(n)
    } else {
      stopifnot(length(w) == n)
    }
    W <- matrix(w, ncol = n, nrow = m, byrow = TRUE)
    Y <- sqrt(1 - rho) * Z + sqrt(rho) * W
  }
  attr(Y, "w") <- w
  Y
}

###########################################################################
# HISTORY:
# 2026-11-25
# o Renaming: camel_case and shorter names
# 2014-04-21
# o Created from 'simulateFactorModelNullsFromSingularValuesAndLoadings'.
###########################################################################

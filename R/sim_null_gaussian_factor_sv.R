#' sim_null_gaussian_factor_sv
#'
#' Simulate null hypotheses according to a factor model
#'
#' @param m Number of tests
#' @param h A vector of \code{k} singular values associated to each factor
#' @param P A \code{m x k} matrix of factor loadings
#' @param rho \code{1-rho} is the standard deviation of the noise
#' @return \item{Y}{a vector of \code{m} simulated observations} \item{W}{a
#' vector of \code{k} factors}
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom Matrix Matrix
#' @importFrom stats rnorm
#' @export
#' @keywords internal
#' @examples
#'
#' m <- 10
#'
#' ## independent
#' sim0 <- sanssouci:::sim_null_gaussian_factor_sv(m)
#' str(sim0)
#' sum(is.na(sanssouci:::sim_null_gaussian_factor_sv(m)))
#' S0 <- sigma_factor_model(m)
#' image(S0)
#'
#' ## equi-correlated
#' rho <- 0.2
#' h <- 1
#' P <- matrix(1, m, length(h))
#' sim1 <- sanssouci:::sim_null_gaussian_factor_sv(m,
#'           h, P, rho=rho)
#' str(sim1)
#' S1 <- sigma_factor_model(m, h, P, rho=rho)
#' image(S1)
#'
#' ## 3-factor model
#' m <- 4*floor(m/4) ## make sure m/4 is an integer
#' rho <- 0.5
#' h <- c(0.5, 0.3, 0.2)*m
#' gamma1 <- rep(1,m)
#' gamma2 <- rep(c(-1, 1), each=m/2)
#' gamma3 <- rep(c(-1, 1), times=m/2)
#' P <- cbind(gamma1, gamma2, gamma3)/sqrt(m)
#'
#' sim3 <- sim_null_gaussian_factor_sv(m,
#'               h, P, rho=rho)
#' str(sim3)
#' S3 <- sigma_factor_model(m, h, P, rho=rho)
#' image(S3)
#' sim3
#'
sim_null_gaussian_factor_sv <- function(m, h = numeric(0),
                                        P = Matrix(0, nrow = m,
                                                   ncol = length(h)),
                                        rho = 0) {
  k <- length(h)
  ## sanity checks
  stopifnot(nrow(P) == m)
  stopifnot(ncol(P) == k)
  if (k > 1) {
    ## is P'P orthonormal ?
    mm <- max(abs(diag(k) - t(P) %*% P))
    if (mm > 1e-10) {
      stop("t(P) %*% P should be orthonormal")
    }
    rm(mm)
  }

  Z <- matrix(rnorm(m), ncol = 1)
  W <- rnorm(k)
  if (rho == 0) {
    Y <- Z
  } else {
    Y <- sqrt(1 - rho) * Z + sqrt(rho) * P %*% Matrix(sqrt(h) * W, ncol = 1)
    #    [code] : if (sum(is.na(Y))) browser()
  }
  Y <- as.vector(Y)
  list(Y = Y, W = W)
}

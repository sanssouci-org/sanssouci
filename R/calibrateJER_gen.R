#' Calibration of joint Family-Wise Error Rate thresholds Generic
#'
#' Calibration of of JER thresholds from one or two-sample tests
#'
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations.
#'   The column names of X should be "0" for the first sample and "1" for the
#'   second sample.
#' @param B A numeric value, the number of permutations to be performed
#' @param alpha Target JER level
#' @param refFamily A character value which can be \describe{
#'
#'   \item{Simes}{The classical family of thresholds introduced by Simes (1986):
#'   \eqn{\alpha*k/m}. This family yields joint FWER control if the test
#'   statistics are positively dependent (PRDS) under H0.}
#'
#'   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k,
#'   \eqn{(t_k)} controls the (marginal) k-FWER.}}
#'
#' @param K For JER control over \code{1:K}, ie joint control of all?
#'   \eqn{k}-FWER, \eqn{k \le K}.
#'   
#' @param test A string indicating the statistical test that should be used. Defaults is 't.test'.
#'   
#' @param parallel If TRUE, the permutations are run in parallel.
#' 
#' @param core If parallel is TRUE, the number of cores to run the permutations in parallel.
#' 
#' @param verbose A boolean value: should extra info be printed?
#' @return A list with elements: \describe{
#'
#'   \item{stat}{A numeric vector of \code{m} test statistics}
#'
#'   \item{thr}{A numeric vector of length \code{K}, such that the estimated
#'   probability that there exists an index \eqn{k} between 1 and \eqn{K} such
#'   that the \eqn{k}-th maximum of the test statistics of is greater than
#'   \eqn{thr[k]}, is less than \eqn{\alpha}}
#'
#'   \item{lambda}{A numeric value, the result of the calibration} }
#'
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#'
#' m <- 543
#' pi0 <- 0.8
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100,
#'                        pi0 = pi0, SNR = 3, prob = 0.5)
#' X <- sim$X
#' cal <- calibrateJER(X, B = 1e3, alpha = 0.2, refFamily="Simes")
#' cal$lambda # > alpha (whp) is rho > 0
#' 
#' # Application 1: confidence envelope
#' #   ie upper confidence bound for the number of false positives 
#' #   among the k most significant items for all k
#' stat <- cal$stat
#' o <- order(stat, decreasing = TRUE)
#' R <- seq_along(stat)
#' Vbar <- curveMaxFP(stat[o], cal$thr)
#' 
#' # True number of false positives among the most significant items:
#' H0 <- which(sim$H == 0)
#' V <- cumsum(o %in% H0)
#' 
#' plot(R, Vbar, t = 's', xlab = "Number of rejections", 
#'                        ylab = "Bound on the number of false positives")
#' lines(R, V, t = 's', col = 2)
#' legend("topleft", c("True number of FP", "Post hoc upper bound"), col = c(2, 1), lty = rep(1, 2))
#' abline(a = 0, b = 1, lty = 2)
#' 
#' # Application 2a: bound on the number of false positives in one or 
#' #    more user-defined selections
#' 
#' sel <- stat[o][c(1:10, 35:40)]
#' maxFP(sel, cal$thr)
#' 
#' sel <- stat[o][c(30:50)]
#' maxFP(sel, cal$thr)
#' 
#' # Application 2b: bound on pi0, the proportion of false positives in H
#' 
#' sel <- stat[o]
#' maxFP(sel, cal$thr)/m
#' pi0
#' 
calibrateJER_gen <- function(X, B, alpha, refFamily = c("Simes", "kFWER"),
                             K = nrow(X), test=NULL, parallel=FALSE,core=NULL, verbose=TRUE) {
  ## sanity checks
  m <- nrow(X);
  refFamily <- match.arg(refFamily)
  
  tests <- testByRandomization_gen(X, B = B, test=test, parallel=parallel, core=core)
  X0 <- tests$T0
  x <- tests$T
  rm(tests)
  
  res <- calibrateJER0(X0, refFamily = refFamily, alpha = alpha, x, kMax = K)
  calib <- list(stat = x, thr = res$thr, lambda = res$lambda) 
  
  return(calib)
}





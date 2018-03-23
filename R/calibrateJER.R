#' Calibration of joint Family-Wise Error Rate thresholds
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
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}.
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
#' sim <- gaussianSamples(m = 123, rho = 0.2, n = 100,
#'                        pi0 = 0.8, SNR = 1, prob = 0.5)
#' X <- sim$X
#' cal <- calibrateJER(X, B = 1e3, alpha = 0.2, refFamily="Simes")
#' cal$lambda # > alpha (whp) is rho > 0

calibrateJER <- function(X, B, alpha, refFamily = c("Simes", "kFWER"),
                            K = nrow(X), verbose=TRUE) {
    ## sanity checks
    m <- nrow(X);
    refFamily <- match.arg(refFamily)

    tests <- testByRandomization(X, B = B)
    
    flavor <- tests$flavor    
    if (flavor == "perm") {
        ## FIXME: why are the below transformation not done internally in the testByRandomization function??
        # p-values
        p0 <- tests$param.p0
        p <- tests$param.p
        # map (back) to the scale of Gaussian test statistics under H0
        X0 <- qnorm(1 - p0) 
        x <- qnorm(1 - p)
    } else {
        X0 <- tests$T0
        x <- tests$T
    }    
    res <- jointFWERControl(X0, refFamily = refFamily, alpha = alpha, x, kMax = K)
    calib <- list(stat = x, thr = res$thr, lambda = res$lambda) 
    
    return(calib)
}



    

#' Calibration of joint Family-Wise Error Rate thresholds
#'
#' Calibration of of JER thresholds using one or two-sample tests
#'
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param categ An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#' @param B A numeric value, the number of permutations to be performed
#' @param alpha Target JER level
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#' @param rowTestFUN A (vectorized) test function. Defaults to
#'   \code{\link{rowWelchTests}}
#' @param refFamily A character value which can be \describe{
#'
#'   \item{Simes}{The classical family of thresholds introduced by Simes (1986):
#'   \eqn{\alpha*k/m}. This family yields joint FWER control if the test
#'   statistics are positively dependent (PRDS) under H0.}
#'
#'    \item{Beta}{A family of thresholds that achieves "balanced" joint error
#'    rate (JER) control under independence}
#'
#'   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k,
#'   \eqn{(t_k)} controls the (marginal) k-FWER.}}
#' @param maxStepsDown Maximum number of steps down to be performed.
#     \code{maxStepsDown=0} corresponds to single step JFWER control. 
#     Defaults to 10.
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}. Defaults to \eqn{m}.
#' @param verbose A boolean value: should extra info be printed?
#' @details See \code{\link{testByRandomization}} for a description of the tests performed for calibration.
#' @return A list with elements: \describe{
#'
#'   \item{p.values}{A numeric vector of \code{m} p-values}
#'   \item{pivStat}{A numeric vector of length \code{B}, the values of the pivotal
#'          statistic whose quantile of order \eqn{alpha} is \eqn{lambda}.}

#'   \item{thr}{A numeric vector of length \code{K}, such that the estimated
#'   probability that there exists an index \eqn{k} between 1 and \eqn{K} such
#'   that the \eqn{k}-th maximum of the test statistics of is greater than
#'   \eqn{thr[k]}, is less than \eqn{\alpha}}
#'
#'   \item{lambda}{A numeric value, the result of the calibration} 
#'   
#'   \item{FP}{A numeric vector of length \code{m}, a 1-alpha confidence envelope on the number of false positives} 
#'   }
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#'
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#'
#' m <- 543
#' pi0 <- 0.8
#' sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
#'                        pi0 = pi0, SNR = 3, prob = 0.5)
#' X <- sim$X
#' categ <- sim$categ
#' alpha <- 0.1
#' cal <- calibrateJER(X, categ, B = 1e2, alpha = alpha, refFamily="Simes")
#' cal$lambda # > alpha (whp) if rho > 0
#' 
#' # Application 1: confidence envelope
#' #   ie upper confidence bound for the number of false positives 
#' #   among the k most significant items for all k
#' env <- cal$conf_env
#' plotConfidenceEnvelope(env, xmax = 200)
#'   
#' ## Compare to Simes (without calibration) and "Oracle" (ie truth from the simulation settings)
#' env_Simes <- confidenceEnvelope(cal$p.values, refFamily = "Simes", param = alpha)
#' env_Oracle <- confidenceEnvelope(cal$p.values, refFamily = "Oracle", param = (sim$H == 0))
#' all_env <- list("Simes + calibration" = env, 
#'                 "Simes"= env_Simes, 
#'                 "Oracle" = env_Oracle)
#' plotConfidenceEnvelope(all_env, xmax = 200)
#' 
#' # Application 2a: bound on the number of false positives in one or 
#' #    more user-defined selections
#' 
#' spval <- sort(cal$p.values)
#' sel <- spval[c(1:10)]
#' maxFP(sel, cal$thr)
#' 
#' sel <- spval[c(1:10, 35:40)]
#' maxFP(sel, cal$thr)
#' 
#' sel <- spval[c(30:50)]
#' maxFP(sel, cal$thr)
#' 
#' # Application 2b: bound on pi0, the proportion of false positives in H
#' 
#' maxFP(spval, cal$thr)/m
#' pi0
#' 
calibrateJER <- function(X, categ, B, alpha, 
                         alternative = c("two.sided", "less", "greater"), 
                         rowTestFUN = rowWelchTests,
                         refFamily = c("Simes", "kFWER", "Beta"),
                         maxStepsDown = 10L,
                         K = nrow(X), verbose=TRUE) {
    alternative <- match.arg(alternative)
    ## sanity checks
    n <- ncol(X)
    m <- nrow(X)
    rownames(X) <- NULL ## make it explicit that rownames are not considered
    stopifnot(length(categ) == n)
    stopifnot(all(categ %in% c(0,1)))
    refFamily <- match.arg(refFamily)
    
    tests <- testByRandomization(X, categ = categ, B = B, alternative = alternative, rowTestFUN = rowTestFUN)
    pval0 <- tests$p0
    pval <- tests$p

    rm(tests)
    res <- calibrateJER0(pval0, refFamily = refFamily, alpha = alpha, 
                         p.values = pval, maxStepsDown = maxStepsDown, kMax = K)
    # fam <- toFamily(refFamily, res$lambda)
    conf_env <- confidenceEnvelope(p.values = pval, refFamily = refFamily, param = res$lambda, K = K)
    proc <- sprintf("%s + calibration", refFamily)
    if (K < m) {
        proc <- sprintf("%s (K = %s)", proc, K)
    }
    conf_env$procedure <- proc
    
    calib <- list(p.values = pval, pivStat = res$pivStat, thr = res$thr, lambda = res$lambda, conf_env = conf_env) 
    return(calib)
}

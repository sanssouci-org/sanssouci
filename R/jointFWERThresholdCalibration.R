##' Joint Family-Wise error rate control
##'
##' Get a family of thresholds that provide joint FWER control
##'
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{ \item{m}{is the number of
##' tested hypotheses} \item{B}{is the number of Monte-Carlo samples}}
##' @param refFamily A character value which can be \describe{ \item{Simes}{The
##' classical family of thresholds introduced by Simes (1986):
##' \eqn{\alpha*k/m}. This family yields joint FWER control if the test
##' statistics are positively dependent (PRDS) under H0.} \item{kFWER}{A family
##' \eqn{(t_k)} calibrated so that for each k, \eqn{(t_k)} controls the
##' (marginal) k-FWER.} }
##' @param alpha Target joint FWER level.
##' @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
##' k[max]}).
##' @param Rcpp If \code{TRUE} (the default), some costly operations (sorting)
##' are performed in C++.
##' @param verbose If \code{TRUE}, print results of intermediate calculations.
##' Defaults to \code{FALSE}.
##'
##' @return A list with elements:\describe{
##' \item{thr}{A numeric vector \code{thr}, such that the estimated
##' probability that there exists an index \eqn{k} between 1 and m
##' such that the k-th maximum of the test statistics of is greater
##' than \eqn{thr[k]}, is less than \eqn{\alpha}.}
##' \item{pivStat}{A numeric vector, the values of the pivotal
##' statistic whose quantile of order \eqn{alpha} is \eqn{lambda}}
##' \item{lambda}{JFWER threshold.}
##' }
##'
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @importFrom stats quantile
##' @export
##' @examples
##'
##' mat <- simulateGaussianNullsFromFactorModel(m=1023, n=1000, flavor="equi-correlated", rho=0.2)
##'
##' alpha <- 0.2
##' resK <- jointFWERControl(mat, refFamily="kFWER", alpha)
##' str(resK)
##'
##' resS <- jointFWERControl(mat, refFamily="Simes", alpha)
##' str(resS)
##'
jointFWERThresholdCalibration <- function(mat, thresholdFamily, pivotalStatFUN, alpha, kMax=nrow(mat), Rcpp=TRUE) {
    ## sanity checks
    m <- nrow(mat);
    stopifnot(length(thresholdFamily(alpha))==m)
    sk <- threshodFamily
    
    ## (single-step) jointFWER threshold calibration
    pivStat <-  pivotalStat(mat, pivotalStatFUN=pivotalStatFUN)
    lambda <- quantile(pivStat, alpha, type=1)
    thr <- sk(lambda)

    res <- list(
        thr=thr, 
        pivStat=pivStat,
        lambda=lambda
    )
}
############################################################################
## HISTORY:
##
## 2016-05-26
## o Created from 'getJointFWERThresholds'.
############################################################################


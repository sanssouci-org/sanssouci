##' Joint Family-Wise error rate threshold calibration
##'
##' Calibrate a family of thresholds that provide joint FWER control
##'
##' This function is mostly intended for internal use.
##'
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{ \item{m}{is the number of
##' tested hypotheses} \item{B}{is the number of Monte-Carlo samples}}
##'
##' @param thresholdFamily A one-parameter reference family as
##' described in BNR. The package implements the following choices:
##' \code{\link{SimesThresholdFamily}} and
##' \code{\link{kFWERThresholdFamily}}.
##'
##' @param pivotalStatFUN A function that takes as input a \eqn{m} x
##' \eqn{B} matrix statistics under H0, and returns the pivotal
##' statistic of BNR. The package implements the following choices:
##' \code{\link{SimesPivotalStatistic}} and \code{\link{SimesPivotalStatistic}}.
##'
##' @param alpha Target joint FWER level.
##' @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
##' k[max]}).
##' @param Rcpp If \code{TRUE} (the default), some costly operations (sorting)
##' are performed in C++.
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
##' @seealso jointFWERControl
##' @examples
##'
##' mat <- simulateGaussianNullsFromFactorModel(m=1023, n=1000, flavor="equi-correlated", rho=0.2)
##' sk <- kFWERThresholdFamily(mat)
##'
##' ## single-step JFWER control
##' alpha <- 0.2
##' res0 <- jointFWERThresholdCalibration(mat,
##'     thresholdFamily=sk, pivotalStatFUN=kFWERPivotalStatistic,
##'     alpha=alpha)
##'
##' ## check coverage
##' prob <- empiricalCoverage(res0$thr, mat);
##' (prob <= alpha)
##'
jointFWERThresholdCalibration <- function(mat, thresholdFamily, pivotalStatFUN, alpha, kMax=nrow(mat), Rcpp=TRUE) {
    ## sanity checks
    m <- nrow(mat)
    stopifnot(kMax<=m)
    stopifnot(length(thresholdFamily(alpha))==m)
    sk <- thresholdFamily

    ## (single-step) jointFWER threshold calibration
    pivStat <-  pivotalStat(mat, kMax=kMax, FUN=pivotalStatFUN)
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


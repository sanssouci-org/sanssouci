##' jointFWERControl
##'
##' get a family of thresholds that provide joint FWER control
##'
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
##' @param \dots Further arguments to be passed to getJointFWERControl.
##' @param verbose If \code{TRUE}, print results of intermediate calculations.
##' Defaults to \code{FALSE}.
##'
##' @return A list with elements:\describe{
##' \item{thr}{A numeric vector \code{thr}, such that the estimated
##' probability that there exists an index \eqn{k} between 1 and m
##' such that the k-th maximum of the test statistics of is greater
##' than \eqn{thr[k]}, is less than \eqn{\alpha}.}
##' \item{prob}{Estimated probability that there exists an index
##' \eqn{k} between 1 and m such that the k-th maximum of the test
##' statistics of is greater than \eqn{thr[k]} (should be in
##' \eqn{[\alpha-tol, alpha]}).}
##' \item{lambda}{JFWER threshold.}
##' \item{sk}{A function such that \code{thr} is identical to
##' \code{sk(lambda)}. If the input argument \code{refFamily} is a
##' function, then \code{sLambda=refFamily}.}
##' \item{pivStat}{A numeric vector, the values of the pivotal
##' statistic whose quantile of order \eqn{alpha} is \eqn{lambda}}
##' \item{Q}{Result of sorting the input score matrix by row and then
##' by columns. Corresponds to matrix 'Q' in Meinshausen (2006).}
##' \item{kMaxH0}{k-max of the test statistics under H0 for \eqn{1
##' \leq k \leq kMax}.}
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
jointFWERControl <- function(mat, refFamily=c("Simes", "kFWER"), alpha, kMax=nrow(mat), Rcpp=TRUE, ..., verbose=FALSE) {
    m <- nrow(mat)
    refFamily <- match.arg(refFamily)
    Q <- NULL
    if (refFamily=="kFWER") {
        sk <- kFWERThresholdFamily(mat, kMax=kMax, Rcpp=Rcpp)
        Q <- attr(sk, 'Q')
    } else if (refFamily=="Simes") {
        sk <- SimesThresholdFamily(m, kMax=kMax)
    }

    pivStat <-  pivotalStat(mat, refFamily=refFamily)
    lambda <- quantile(pivStat, alpha, type=1)
    thr <- sk(lambda)

    kmaxH0 <- partialColSortDesc(mat, nrow(mat));
    prob <- coverage(thr, kmaxH0);

    res <- list(
        thr=thr,
        prob=prob,
        lambda=lambda,
        sk=sk,
        pivStat=pivStat,
        Q=Q,
        kmaxH0=kmaxH0
    )
}
############################################################################
## HISTORY:
##
## 2016-05-26
## o Created from 'getJointFWERThresholds'.
############################################################################


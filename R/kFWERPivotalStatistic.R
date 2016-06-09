##' Pivotal statistic for the balanced threshold family of BNR
##'
##' @param kmaxH0 A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
##' samples of \code{c} test statistics under the null hypothesis,
##' with each sample (column) sorted decreasingly.
##' @param m An integer value, the total number of hypotheses
##' tested. \code{m} should not be less than \code{c}.
##' @return the pivotal statistic \eqn{s_k} for the balanced threshold family introduced by BNR
##'
##' @export
##'
kFWERPivotalStatistic <- function(kmaxH0, m) {
    B <- ncol(kmaxH0)
    rks <- rowRanks(kmaxH0, ties.method="max") ## 'max' is the default in 'matrixStats'
    colMins(rks)/B
}

##' Pivotal statistic for the balanced threshold family of BNR
##'
##' @param kmaxH0 A \eqn{m} x \eqn{B} matrix of \code{B} Monte-Carlo samples of \code{m} test
##' statistics under the null hypothesis, with each sample (column) sorted decreasingly.
##'
##' @return the pivotal statistic \eqn{s_k} for the balanced threshold family introduced by BNR
##'
##' @export
##'
kFWERPivotalStatistic <- function(kmaxH0) {
    B <- ncol(kmaxH0)
    rks <- rowRanks(kmaxH0, ties.method="max") ## 'max' is the default in 'matrixStats'
    colMins(rks)/B
}

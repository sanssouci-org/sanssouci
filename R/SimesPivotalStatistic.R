##' Simes' pivotal statistic
##'
##' The pivotal statistic for Simes' threshold family
##'
##' @param kmaxH0 A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
##' samples of \code{c} test statistics under the null hypothesis,
##' with each sample (column) sorted decreasingly.
##' @param m An integer value, the total number of hypotheses
##' tested. \code{m} should not be less than \code{c}.
##' @return the pivotal statistic \eqn{s_k} of BNR for Simes'
##' threshold family
##'
##' @details The inverse function of \eqn{s_k} is: \eqn{s_k^{-1}:u
##' \mapsto (m/k)*(1-pnorm(u))}.
##'
##' @export
##'
SimesPivotalStatistic <- function(kmaxH0, m) {
    ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
    c <- nrow(kmaxH0)
    stopifnot(c<=m)
    pval <- 1-pnorm(kmaxH0)
    skInv <- sweep(pval, 1, m/1:c, "*")
    colMins(skInv)
}

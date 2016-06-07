##' Simes' pivotal statistic
##'
##' The pivotal statistic for Simes' threshold family
##'
##' @param kmaxH0 A \eqn{m} x \eqn{B} matrix of \code{B} Monte-Carlo samples of \code{m} test
##' statistics under the null hypothesis, with each sample (column) sorted decreasingly.
##'
##' @return the pivotal statistic \eqn{s_k} of BNR for Simes' threshold family
##'
##' @details The inverse function of \eqn{s_k} is: \eqn{s_k^{-1}:u
##' \mapsto (m/k)*(1-pnorm(u))}.
##'
##' @export
##'
SimesPivotalStatistic <- function(kmaxH0) {
    ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
    m <- nrow(kmaxH0)
    pval <- 1-pnorm(kmaxH0)
    skInv <- sweep(pval, 1, m/1:m, "*")
    colMins(skInv)
}

# Simes' pivotal statistic
# 
# The pivotal statistic for Simes' threshold family
# 
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} samples of \code{c} test
#   statistics under the null hypothesis.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq 
#   k[max]}).
# @param m An integer value, the total number of hypotheses tested. \code{m}
#   should not be less than \code{c}.
# @return the pivotal statistic \eqn{s_k} of BNR for Simes' threshold family
#   
# @details The inverse function of \eqn{s_k} is: \eqn{s_k^{-1}:u \mapsto
#   (m/k)*(1-pnorm(u))}.
#   
#' @importFrom matrixStats colMins
#' @importFrom stats pnorm
#   
SimesPivotalStatistic <- function(mat, kMax, m) {
    stopifnot(kMax <= m)
    ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
    B <- ncol(mat)
    c <- min(kMax, nrow(mat))  # K \vee |C| in the BNR paper
    
    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    kmaxH0 <- partialColSortDesc(mat, c);
    pval <- 1 - pnorm(kmaxH0)
    skInv <- sweep(pval, MARGIN = 1, STATS = m/1:c, FUN = "*")
    colMins(skInv)
}

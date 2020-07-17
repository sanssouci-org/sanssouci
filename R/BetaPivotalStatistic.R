# Beta pivotal statistic
# 
# The pivotal statistic for the Beta threshold family
# 
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} samples of \code{c} test
#   statistics under the null hypothesis.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq 
#   k[max]}).
# @param m An integer value, the total number of hypotheses tested. \code{m}
#   should not be less than \code{c}.
# @return the pivotal statistic \eqn{s_k} of BNR for the Beta threshold family
#   
# @details The inverse function of \eqn{s_k} is: \eqn{s_k^{-1}:u \mapsto
#   pbeta(1-pnorm(u), k, m + 1 - k))}.
#   
#' @importFrom matrixStats colMins
#' @importFrom stats pbeta pnorm
#   
BetaPivotalStatistic <- function(mat, kMax, m) {
    stopifnot(kMax <= m)
    ## for Beta, s_k^{-1}(u) = pbeta(1-pnorm(u), k, m + 1 - k))
    c <- min(kMax, nrow(mat))  # K \vee |C| in the BNR paper
    
    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    kmaxH0 <- partialColSortDesc(mat, c);
    pval <- pnorm(kmaxH0, lower.tail = FALSE)
    skInv <- matrix(nrow = nrow(pval), ncol = ncol(pval))
    for (kk in 1:kMax) {
        skInv[kk, ] <- pbeta(pval[kk, ], kk, m + 1 - kk)
    }
    colMins(skInv)
}

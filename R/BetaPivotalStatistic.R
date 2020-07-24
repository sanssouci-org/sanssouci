# Beta pivotal statistic
# 
# The pivotal statistic for the Beta threshold family
# 
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} samples of \code{c} 
#   p-values under the null hypothesis.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq 
#   k[max]}).
# @param m An integer value, the total number of hypotheses tested. \code{m}
#   should not be less than \code{c}.
# @return the pivotal statistic of BNR for the Beta threshold family
#   
# @details The inverse function of \eqn{t_k} is: \eqn{t_k^{-1}:u \mapsto
#   pbeta(u, k, m + 1 - k))}.
#   
#' @importFrom matrixStats colMins
#' @importFrom stats pbeta
#   
BetaPivotalStatistic <- function(mat, kMax, m) {
    stopifnot(kMax <= m)
    ## for Beta, t_k^{-1}(u) = pbeta(u, k, m + 1 - k))
    c <- min(kMax, nrow(mat))  # K \vee |C| in the BNR paper
    
    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    pval <- -partialColSortDesc(-mat, c);
    tkInv <- matrix(nrow = nrow(pval), ncol = ncol(pval))
    for (kk in 1:kMax) {
        tkInv[kk, ] <- pbeta(pval[kk, ], kk, m + 1 - kk)
    }
    colMins(tkInv)
}

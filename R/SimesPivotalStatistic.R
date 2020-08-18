# Simes' pivotal statistic
# 
# The pivotal statistic for Simes' threshold family
# 
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} samples of \code{c} p-values
#  under the null hypothesis.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq 
#   k[max]}).
# @param m An integer value, the total number of hypotheses tested. \code{m}
#   should not be less than \code{c}.
# @return the pivotal statistic of BNR for Simes' threshold family
#   
# @details The inverse function of \eqn{t_k} is: \eqn{t_k^{-1}:u \mapsto
#   (m/k)*u}.
#   
#' @importFrom matrixStats colMins
#   
SimesPivotalStatistic <- function(mat, kMax, m) {
    stopifnot(kMax <= m)
    ## for Simes, t_k^{-1}(u) = (m/k)*u
    c <- min(kMax, nrow(mat))  # K \vee |C| in the BNR paper
    
    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    pval <- -partialColSortDesc(-mat, c);
    tkInv <- sweep(pval, MARGIN = 1, STATS = m/1:c, FUN = "*")
    colMins(tkInv)
}

##' Pivotal statistic for the balanced threshold family of BNR
##'
##' @param kmaxH0 A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
##' samples of \code{c} test statistics under the null hypothesis,
##' with each sample (column) sorted decreasingly.
##' @param m An integer value, the total number of hypotheses
##' tested. \code{m} should not be less than \code{c}.
##' @return the pivotal statistic \eqn{s_k} for the balanced threshold family introduced by BNR
##'
kFWERPivotalStatistic <- function(mat, kMax=nrow(mat), C=1:nrow(mat)) {
    c <- length(C)
    c <- min(kMax, c)  # K \vee |C| in th BNR paper

##    kmaxH0 <- partialColSortDesc(mat, k=kMax);
    kmaxH0 <- partialColSortDesc(mat, c);  ## no need to go further than c!
    kmaxH0C <- partialColSortDesc(mat[C, ], c);
    minPseudoRanks(kmaxH0, kmaxH0C)
}

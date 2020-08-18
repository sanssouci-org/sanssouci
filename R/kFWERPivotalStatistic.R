# Pivotal statistic for the balanced threshold family of BNR
#
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
# samples of \code{c} p-values under the null hypothesis,
# with each sample (column) sorted decreasingly.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq
# k[max]}).
# @param C A subset of \code{1:c} on which calibration is performed.
# @return the pivotal statistic \eqn{s_k} for the balanced threshold family introduced by BNR
#
kFWERPivotalStatistic <- function(mat, kMax=nrow(mat), C=1:nrow(mat)) {
    c <- length(C)
    c <- min(kMax, c)  # K \vee |C| in th BNR paper: no need to go further than c!

##    kmaxH0 <- partialColSortDesc(mat, k=kMax);
    kmaxH0 <- -partialColSortDesc(-mat, c);
    if (identical(C, 1:nrow(mat))) { ## avoid re-calculating the same statistics twice...
        kmaxH0C <- kmaxH0
    } else {
        kmaxH0C <- -partialColSortDesc(-mat[C, ], c);
    }
    minPseudoRanks(-kmaxH0, -kmaxH0C)
}

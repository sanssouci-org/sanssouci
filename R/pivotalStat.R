##' Pivotal statistic of BNR
##'
##' Calculate the pivotal statistic of BNR
##'
##' @param mat A \eqn{c} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{
##'   \item{c}{is the number of null hypotheses tested (possibly less than \code{m})}
##'   \item{B}{is the number of Monte-Carlo samples}}
##' @param FUN A function implementing the calculation of pivotal
##' statistic for a specific reference family. Can either by
##' \code{kFWERPivotalStatistic} or \code{SimesPivotalStatistic}.
##' @param m An integer value, the total number of hypotheses
##' tested. \code{m} should not be less than \code{c}.
##' @param kMax A scalar value between \code{1} and \code{m} such that
##' simultaneous control of (\eqn{k}-FWER for all \eqn{k \le k[max]}) is
##' targeted.
##' @param Rcpp If \code{TRUE} (the default), some costly operations
##' (sorting) are performed in C++.
##' @return A numeric vector, the values of the pivotal statistic
##' whose quantile of order \eqn{alpha} is the \eqn{lambda}-adjustment
##' factor of BNR
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @importFrom stats pnorm
##' @importFrom matrixStats colMins rowRanks
##' @examples
##'
##' mat <- simulateGaussianNullsFromFactorModel(m=1023, n=1000, flavor="equi-correlated", rho=0.2)
##' piv <- pivotalStat(mat, FUN=kFWERPivotalStatistic)
##' str(piv)
##'
##' ## check JFWER control using these statistics:
##' alpha <- 0.2
##' lambda <- quantile(piv, alpha, type=1)
##'
##' sk <- kFWERThresholdFamily(mat)
##' thr <- sk(lambda)
##' prob <- empiricalCoverage(thr, mat);
##' (prob <= alpha)
##'
##'
pivotalStat <- function(mat,
                        FUN,
                        m=nrow(mat),
                        kMax=nrow(mat),
                        Rcpp=TRUE){
    B <- ncol(mat)
    c <- min(kMax, nrow(mat))  # K \vee |C| in th BNR paper
    stopifnot(kMax<=m)

    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    if (Rcpp) {
        kmaxH0 <- partialColSortDesc(mat, c);
    } else {
        kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)  ## k-max of the test statistics under H0
        kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]          ## truncate to [1,kMax]
    }

    pivotalStat <- FUN(kmaxH0, m)
    return(pivotalStat)
}


############################################################################
## HISTORY:
##
## 2016-05-24
## o Created from getJointFWERThresholds.
############################################################################


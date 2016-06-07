##' Pivotal statistic of BNR
##'
##' Calculate the pivotal statistic of BNR
##'
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{ \item{m}{is the number of
##' null hypotheses tested} \item{B}{is the number of Monte-Carlo samples}}
##' @param FUN A character value which can be \describe{ \item{Simes}{The
##' classical family of thresholds introduced by Simes (1986):
##' \eqn{\alpha*k/m}. This family yields joint FWER control if the test
##' statistics are positively dependent (PRDS) under H0.} \item{kFWER}{A family
##' \eqn{(t_k)} calibrated so that for each k, \eqn{(t_k)} controls the
##' (marginal) k-FWER.} }
##' @param kMax A scalar value between \code{1} and \code{m} such that
##' simultaneous control of (\eqn{k}-FWER for all \eqn{k \le k[max]}) is
##' targeted.
##' @param Rcpp If \code{TRUE}, some costly operations (sorting) are performed
##' in C++.
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
                        kMax=nrow(mat),
                        Rcpp=FALSE){
    m <- nrow(mat)
    B <- ncol(mat)

    ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
    if (Rcpp) {
        kmaxH0 <- partialColSortDesc(mat, kMax);
    } else {
        kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)        ## k-max of the test statistics under H0:
        kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]                ## truncate to [1,kMax]
    }

    pivotalStat <- FUN(kmaxH0)
    return(pivotalStat)
}


############################################################################
## HISTORY:
##
## 2016-05-24
## o Created from getJointFWERThresholds.
############################################################################


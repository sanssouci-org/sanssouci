##' # calculate the pivotal statistic of BNR.
##'
##'
##'
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{ \item{m}{is the number of
##' null hypotheses tested} \item{B}{is the number of Monte-Carlo samples}}
##' @param refFamily A character value which can be \describe{ \item{Simes}{The
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
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @importFrom stats pnorm
##' @importFrom matrixStats colMins rowRanks
##' @examples
##'
##' m <- 1023
##' B <- 1e3
##'
##' flavor <- c("independent", "equi-correlated", "3-factor model")[2]
##' rho <- 0.2
##' mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
##' piv <- pivotalStat(mat, refFamily="kFWER")
##' pivS <- pivotalStat(mat, refFamily="Simes")
##'
##' ## check JFWER control using these statistics:
##' alpha <- 0.2
##' lambda <- quantile(piv, alpha, type=1)
##' Q <- bisort(mat, Rcpp=TRUE)
##' sk <- function(alpha) thresholdFamily(alpha, Q)
##' thr <- sk(lambda)
##' kmaxH0 <- partialColSortDesc(mat, nrow(mat));
##' prob <- coverage(thr, kmaxH0);
##'
##'
pivotalStat <- structure(function( ## calculate the pivotal statistic of BNR.
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
### statistics under the null hypothesis. \describe{ \item{m}{is the
### number of null hypotheses tested} \item{B}{is the number of
### Monte-Carlo samples}}
    refFamily=c("Simes", "kFWER"),
### A character value which can be \describe{ \item{Simes}{The
###   classical family of thresholds introduced by Simes (1986):
###   \eqn{\alpha*k/m}. This family yields joint FWER control if the
###   test statistics are positively dependent (PRDS) under H0.}
###   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k,
###   \eqn{(t_k)} controls the (marginal) k-FWER.}
### }
    kMax=nrow(mat),
### A scalar value between \code{1} and \code{m} such that
### simultaneous control of (\eqn{k}-FWER for all \eqn{k \le k[max]})
### is targeted.
    Rcpp=FALSE)
### If \code{TRUE}, some costly operations (sorting) are performed in C++.
    {
        m <- nrow(mat)
        B <- ncol(mat)
        refFamily <- match.arg(refFamily)

        ## get matrix 'M' of BNR by (partial) sorting of hypotheses for each sample
        if (Rcpp) {
            kmaxH0 <- partialColSortDesc(mat, kMax);
        } else {
            kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)        ## k-max of the test statistics under H0:
            kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]                ## truncate to [1,kMax]
        }

        ## get pivotal statistic
        if (refFamily=="Simes") {
            ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
            pval <- 1-pnorm(kmaxH0)
            skInv <- sweep(pval, 1, m/1:m, "*")
            pivotalStat <- colMins(skInv)
        } else if (refFamily=="kFWER") {
            rks <- rowRanks(kmaxH0, ties.method="max") ## 'max' is the default in 'matrixStats'
            pivotalStat <- colMins(rks)/B
        }
        return(pivotalStat) ##<< A numeric vector, the values of the pivotal statistic whose quantile of order \eqn{alpha} is the \eqn{lambda}-adjustment factor of BNR
    }, ex=function(){
        m <- 1023
        B <- 1e3

        flavor <- c("independent", "equi-correlated", "3-factor model")[2]
        rho <- 0.2
        mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
        piv <- pivotalStat(mat, refFamily="kFWER")
        pivS <- pivotalStat(mat, refFamily="Simes")

        ## check JFWER control using these statistics:
        alpha <- 0.2
        lambda <- quantile(piv, alpha, type=1)
        Q <- bisort(mat, Rcpp=TRUE)
        sk <- function(alpha) referenceFamily(alpha, Q)
        thr <- sk(lambda)
        kmaxH0 <- partialColSortDesc(mat, nrow(mat));
        prob <- coverage(thr, kmaxH0);

    })


############################################################################
## HISTORY:
##
## 2016-05-24
## o Created from getJointFWERThresholds.
############################################################################


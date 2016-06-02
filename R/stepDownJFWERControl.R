##' stepDownJointFWERControl
##'
##' get a family of thresholds that provide joint FWER control
##'
##' At 'step 0', a one-parameter family \code{sk} ensuring JFWER control at
##' level \code{alpha} is inferred using \code{getJointFWERThresholds}. This
##' family is then kept fixed throughout the step-down process. What changes
##' during the steps down is the over-estimation of \eqn{H_0}: starting with
##' \eqn{H} at step 0, a first lower bound is derived from \code{sk(lambda)},
##' which allows re-calibration of \eqn{lambda} by applying
##' \code{getJointFWERThresholds} to the subset of not-yet-rejected hypotheses,
##' etc.  Note that this family does not change during the step-down process!
##' Only the value of \code{lambda} may change
##'
##' @param stat A vector of \eqn{m} test statistics.
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
##' statistics under the null hypothesis. \describe{ \item{m}{is the number of
##' tested hypotheses} \item{B}{is the number of Monte-Carlo samples}}
##' @param refFamily A character value or a function that returns a vector of
##' \code{m} thresholds (see Details).
##' @param alpha Target joint FWER level.
##' @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
##' k[max]}).
##' @param maxSteps Maximum number of steps down to be performed.
##' \code{maxSteps=1} corresponds to single step JFWER control.
##' @param Rcpp If \code{TRUE}, some costly operations (sorting) are performed
##' in C++.
##' @param \dots Further arguments to be passed to
##' \code{\link{jointFWERControl}}.
##' @param verbose If \code{TRUE}, print results of intermediate calculations.
##' Defaults to \code{FALSE}.
##' @return List with elements: \describe{
##' \item{thr}{A numeric vector \code{thr} of length \eqn{m}, such
##' that the estimated probability that there exists an index \eqn{k}
##' between 1 and m such that the \eqn{k}-th maximum of the test
##' statistics of is greater than \eqn{thr[k]}, is less than
##' \eqn{\alpha}.}xs
##' \item{thrMat}{A matrix of size \eqn{m} x the number of steps down,
##' containing successive JFWER threshold estimates.}
##' \item{pivMat}{A
##' matrix of size \eqn{m} x the number of steps down, containing
##' successive pivotal test statistics}
##' \item{sk}{threshold family}
##' }
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @examples
##'
##' set.seed(0xBEEF)
##' sim <- simulateMein2006(m=1e2, rho=0.2, n=300, pi0=0.9, SNR=2)
##' X <- sim$X
##' y <- sim$y
##' H0 <- which(sim$H==0)
##' H1 <- which(sim$H==1)
##' m0 <- length(H0)
##' m1 <- length(H1)
##'
##' ## Test statistics
##' w <- wilcoxStat(X, y, B=ncol(X))
##' scoreMat <- w$stat0Mat
##' stat <- w$stat
##'
##' pch <- 20
##' ## show test statistics
##' plot(stat, col=rep(c(1, 2), times=c(m0, m1)), main="Test statistics", pch=pch)
##' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
##'
##' alpha <- 0.1
##' resSD <- stepDownJointFWERControl(stat, scoreMat, refFamily="kFWER", alpha=alpha, verbose=TRUE)
##'
##' ## confidence envelopes
##' thr <- resSD$thr
##' o <- order(stat, decreasing=TRUE)
##' statO <- stat[o]
##'
##' Vbar <- upperBoundFP(statO, thr)
##'
##' thrMat <- resSD$thrMat
##' bounds <- apply(thrMat, 2, function(thr) upperBoundFP(statO, thr, flavor="Mein2006"))
##' ## True V (number of false discoveries among first rejections)
##' V <- cumsum(o %in% H0)
##'
##'
stepDownJointFWERControl <- function(
    stat,
    mat,
    refFamily=c("Simes", "kFWER"),
    alpha,
    kMax=nrow(mat),
    maxSteps=Inf,
    Rcpp=FALSE,
    ...,
    verbose=FALSE)
    {
        m <- nrow(mat)
        stopifnot(length(stat)==m)
        refFamily <- match.arg(refFamily)

        ## get threshold family
        Q <- NULL
        if (refFamily=="kFWER") {
            if (is.null(mat)) {
                stop("Argument 'mat' should be provided for family 'kFWER'")
            }
            sk <- kFWERThresholdFamily(mat, kMax=kMax, Rcpp=Rcpp)
            Q <- attr(sk, 'Q')
            attr(sk, 'Q') <- NULL
        } else if (refFamily=="Simes") {
            if (!is.null(mat)) {
                warning("Unused argument 'mat' for family 'Simes'")
            }
            sk <- SimesThresholdFamily(m, kMax=kMax)
        }
        ## single-step thresholds
        pivStat <-  pivotalStat(mat, refFamily=refFamily)
        lambda <- quantile(pivStat, alpha, type=1)
        thr <- sk(lambda)

        ## for storing results
        thrMat <- thr
        pivMat <- pivStat
        converged <- FALSE

        ## step 1
        step <- 1
        thr1 <- thr[1]   ## (1-)FWER threshold
        R1 <- which(stat>=thr1)

        while (!converged && step<maxSteps) {
            step <- step+1

            if (verbose) {
                print(paste("Step:", step))
                print(paste("|R1|=", length(R1)))
            }

            if (length(R1)) {
                ## updated score matrix given R1
                if (verbose) {
                    print(R1)
                }
                mat1 <- mat[-R1, ]
                sk1 <- function(alpha) sk(alpha)[-R1]
            } else {
                mat1 <- mat
                sk1 <- sk
            }

            ## joint FWER control through lambda-adjustment, *holding sk fixed*
            pivStat <-  pivotalStat(mat1, refFamily=refFamily)
            lambda <- quantile(pivStat, alpha, type=1)
            thr <- sk1(lambda)
            thr <- c(thr, rep(-Inf, length(R1)))
            thrMat <- cbind(thrMat, thr)
            pivMat <- cbind(pivMat, pivStat)

            thr1 <- thr[1]   ## (1-)FWER threshold
            R1new <- which(stat>=thr1)

            ## convergence reached?
            converged <- identical(R1new, R1)

            ## update R1
            R1 <- R1new
        }

        ## TODO: also return step at which hyp j is rejected!

        list(thr=thr,
             thrMat=thrMat,
             pivMat=pivMat,
             sk=sk)
    }

############################################################################
## HISTORY:
##
## 2016-02-05
## o de-Implemented the "Oracle" version of step-down control, as it
## can be obtained direclty by tweaking the 'stat' parameter (see
## example).
##
## 2016-01-07
## o BUG FIX: the one-parameter family should not be updated at each step down!
## o Example fixed accordingly.
## o Implemented the "Oracle" version of step-down control.
##
## 2014-05-09
## o Created from 'getJointFWERThresholds'.
############################################################################


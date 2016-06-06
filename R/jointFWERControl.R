##' Calibration of joint Family-Wise error rate thresholds
##'
##' Calibration of a family of thresholds that provide joint FWER control
##'
##' @param stat A vector of \eqn{m} test statistics.
##' @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of
##'     test statistics under the null hypothesis. \describe{
##'     \item{m}{is the number of tested hypotheses} \item{B}{is the
##'     number of Monte-Carlo samples}}
##' @param referenceFamily A character value, the reference family for
##'     calibration (see details).  \code{m} thresholds (see Details).
##' @param alpha Target joint FWER level.
##' @param kMax For simultaneous control of (\eqn{k}-FWER for all
##'     \eqn{k \le k[max]}).
##' @param maxSteps Maximum number of steps down to be performed.
##'     \code{maxSteps=1} corresponds to single step JFWER control.
##' @param Rcpp If \code{TRUE}, some costly operations (sorting) are
##'     performed in C++.
##' @param \dots Further arguments to be passed to
##'     \code{\link{jointFWERControl}}.
##' @return A list with elements: \describe{

##'     \item{thr}{A numeric vector \code{thr}, such that the
##'     estimated probability that there exists an index \eqn{k}
##'     between 1 and m such that the k-th maximum of the test
##'     statistics of is greater than \eqn{thr[k]}, is less than
##'     \eqn{\alpha}.}

##'     \item{pivStat}{A numeric vector, the values of the pivotal
##'     statistic whose quantile of order \eqn{alpha} is \eqn{lambda}.}

##'     \item{lambda}{JFWER threshold.}  

##'     \item{Vbar}{An upper bound on the number of false discoveries,
##'     as calculated by \code{upperBoundFP(stat, thr)}.}

##'     \item{steps}{a list with elements named 'thr', 'pivStat' and
##'     'lambda' giving the sequence of corresponding vectors/values
##'     along the steps down.}}
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
##' res <- jointFWERControl(stat, scoreMat, refFamily="kFWER", alpha=alpha)
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

jointFWERControl <- function(stat=NULL,
                             mat=NULL,
                             refFamily=c("Simes", "kFWER"),
                             kMax=nrow(mat),
                             maxStepsDown=100,
                             alpha,
                             Rcpp=TRUE) {
    ## This function is the main workhorse of the package.

    ## sanity checks
    m <- nrow(mat);
    if (is.null(stat)) {
        stopifnot(refFamily=="Simes")
    } else {
        stopifnot(length(stat)==m)
    }
    stopifnot(kMax<=m)

    if (refFamily=="Simes") {
        sk <- SimesThresholdFamily(m, kMax=m)
        pivStatFUN <- SimesPivotalStatistic
    } else if (refFamily=="kFWER") {
        sk <- kFWERThresholdFamily(mat, kMax=kMax, Rcpp=Rcpp)
        pivStatFUN <- kFWERPivotalStatistic
    }

    ## single-step JFWER control
    res0 <- jointFWERThresholdCalibration(mat, thresholdFamily=sk, pivotalStatFUN=pivotalStatFUN,
                             alpha=alpha, kMax=kMax, Rcpp=Rcpp)
    lambda <- res0$lambda
    thr <- sk(lambda)

    ## storing results
    thrMat <- matrix(thr, ncol=1)
    pivMat <- matrix(pivStat, ncol=1)
    lambdas <- lambda
    converged <- FALSE

    ## step 0
    step <- 0
    thr1 <- thr[1]   ## (1-)FWER threshold
    R1 <- which(stat>=thr1)

    while (!converged && step<maxSteps) {
        step <- step+1

        ## backup
        lambda0 <- lambda
        thr0 <- thr
        R10 <- R1

        if (length(R1)) {
            matL1 <- matL[-R1, ]
            sk1 <- function(alpha) sk(alpha)[-R1]
        } else {
            matL1 <- matL
            sk1 <- sk
        }

        ## joint FWER control through lambda-adjustment, *holding sk fixed*
        res <- jointFWERControl(mat1, sk1, pivotalStatFUN, alpha, kMax=kMax, Rcpp=Rcpp)
        thr <- res$thr
        pivStat <-  res$pivStat
        lambda <-res$lambda1
        thr <- c(thr, rep(-Inf, length(R1)))

        thr1 <- thr[1]   ## (1-)FWER threshold
        R1 <- which(stat>=thr1)

        ## convergence reached?
        noNewRejection <- all(R1 %in% R10)
        ## In rare situations R1 is strictly included in R10,
        ## we need to declare convergence in such cases as well
        ## (and keep the largest rejection set!)
        if (noNewRejection) {
            if (!identical(R1, R10)) {
                ## not a 'TRUE' convergence: override the last step down!
                thr <- thr0
                lambda <- lambda0
            }
            converged <- TRUE ## stop the step-down process
        }

        thrMat <- cbind(thrMat, thr)
        pivMat <- cbind(pivMat, pivStat)
        lambdas <- c(lambdas, lambda)
    }
    if (step==maxSteps) {
        warning("Maximal number of steps down reached without reaching convergence")
    }

    ## upper bound on the number of false positives among first 'natural' rejections
    o <- order(stat, decreasing=TRUE)
    Vbar <- upperBoundFP(stat[o], thr)
    
    stepsDown <- list(
        thr=thrMat,
        pivStat=pivMat,
        lambda=lambdas)
    res <- list(thr=thr,
                pivStat=pivStat,
                lambda=lambda, 
                Vbar=Vbar,
                stepsDown=stepsDown)
    ## TODO: also return step at which hyp j is rejected as in Romano-Wolf's programs?
    ## TODO: also return 'Sbar' (aka 'lowerBoundTP')
    return(res);    
}

                                         

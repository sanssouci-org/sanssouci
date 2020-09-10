# Calibration of joint Family-Wise error rate thresholds
#
# Calibration of of JER thresholds from null test statistics
#
# @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of
#     p-values under the null hypothesis. \describe{
#     \item{m}{is the number of tested hypotheses} \item{B}{is the
#     number of Monte-Carlo samples}}
# @param refFamily A character value which can be \describe{
# \item{Simes}{The classical family of thresholds introduced by
# Simes (1986): \eqn{\alpha*k/m}. This family yields joint FWER
# control if the null hypotheses are positively dependent (PRDS)
# under H0.}
# \item{Beta}{A family of thresholds that achieves "balanced" joint 
#  error rate (JER) control under independence}
# \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k,
# \eqn{(t_k)} controls the (marginal) k-FWER.}}
# @param alpha Target joint FWER level.
# @param p.values A vector of \eqn{m} p.values. Not used for
# single step control, and mandatory for step-down JFWER control. If
# not provided, single step control is performed.
# @param maxStepsDown Maximum number of steps down to be performed.
#     \code{maxStepsDown=0} corresponds to single step JFWER control.
#     Defaults to 10.
# @param kMax For simultaneous control of (\eqn{k}-FWER for all
#     \eqn{k \le k[max]}).
# @param Rcpp If \code{TRUE}, some costly operations (sorting) are
#     performed in C++.
# @param verbose A boolean value: should extra info be printed?
# @return A list with elements: \describe{
#     \item{thr}{A numeric vector of length \code{m}, such that the
#     estimated probability that there exists an index \eqn{k} between 1
#     and m such that the k-th minimum of the p-values less than \eqn{thr[k]}, 
#     is less than \eqn{\alpha}.}
#     \item{pivStat}{A numeric vector of length \code{m}, the values of the pivotal
#     statistic whose quantile of order \eqn{alpha} is \eqn{lambda}.}
#     \item{lambda}{JFWER threshold.}
#     \item{Vbar}{An upper bound on the number of false discoveries,
#     as calculated by \code{curveMaxFP(stat, thr)}.}
#     \item{steps}{a list with elements named 'thr', 'pivStat' and
#     'lambda' giving the sequence of corresponding vectors/values
#     along the steps down.}}
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
# @examples
#
# set.seed(0xBEEF)
# m <- 123
# rho <- 0.2
# n <- 100
# pi0 <- 0.5
# B <- 1e3
# 
# sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = 3)
# X0 <- sim$X0
# x <- sim$x
# 
# ## Test statistics
# pch <- 20
# plot(x, col=1+sim$H, main="Test statistics", pch=pch)
# legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
# 
# alpha <- 0.1
# res <- calibrateJER0(pnorm(X0, lower.tail = FALSE), refFamily="Simes", alpha=alpha, p.values=pnorm(x, lower.tail = FALSE))
# Vbar <- res$Vbar
# ttl <- paste("Upper bound on the number of false positives", "among first k hypotheses", sep="\n")
# plot(Vbar, main=ttl, xlab="k", ylab=expression(bar(V(k))))


calibrateJER0 <- function(mat,
                          refFamily = c("Simes", "kFWER", "Beta"),
                          alpha,
                          p.values = NULL,
                          maxStepsDown = 10L,
                          kMax = nrow(mat),
                          Rcpp = TRUE,
                          verbose = FALSE) {
    ## This function is the main internal workhorse of the package.
    
    ## sanity checks
    m <- nrow(mat);
    refFamily <- match.arg(refFamily)
    if (is.null(p.values)) {
        ## force single step control
        if (verbose && (maxStepsDown > 0)) {
            print("Arguement 'p.values' not provided: cannot perform step-down control")
        }
        maxStepsDown <- 0
    } else {
        stopifnot(length(p.values)==m)
    }
    kMax <- min(kMax, m)  ## 'm' can be greater than 'kMax' throughout the step-down process
    if (verbose) {
        proc <- ifelse(maxStepsDown == 0, "Single step", "Step down")
        msg <- sprintf("Joint Family-Wise Error Rate control: %s procedure based on %s family", proc, refFamily)
        print(msg)
    }
    if (refFamily == "Simes") {
        tk <- SimesThresholdFamily(m, kMax = kMax)
        pivStatFUN <- function(mat, kMax, C) {
            SimesPivotalStatistic(mat[C, ], kMax, nrow(mat))
        }
    } else if (refFamily == "Beta") {
        tk <- BetaThresholdFamily(m, kMax = kMax)
        pivStatFUN <- function(mat, kMax, C) {
            BetaPivotalStatistic(mat[C, ], kMax, nrow(mat))
        }
    } else if (refFamily == "kFWER") {
        tk <- kFWERThresholdFamily(mat, kMax = kMax, Rcpp = Rcpp)
        pivStatFUN <- function(mat, kMax, C) {
            kFWERPivotalStatistic(mat, kMax, C)
        }
    }
    ## (single-step) joint FWER control
    ## pivStat <-  pivotalStat(mat, m=m, kMax=kMax, FUN=pivStatFUN)
    pivStat <-  pivStatFUN(mat, kMax = kMax, 1:m)
    lambda <- stats::quantile(pivStat, alpha, type = 1)
    thr <- tk(lambda)
    
    ## storing results
    thrMat <- matrix(thr, ncol = 1)
    pivMat <- matrix(pivStat, ncol = 1)
    lambdas <- lambda
    
    ## step 0
    step <- 0
    thr1 <- thr[1]   ## (1-)FWER threshold
    if (is.null(p.values)) {
        R1 <- integer(0)
    } else {
        R1 <- which(p.values <= thr1)
    }
    
    
    ## force 'convergence' if nb of "FWER rejections" is 0 (nothing to
    ## gain) or m (nothing left to be rejected)
    converged <- (length(R1) == 0L)  | (length(R1) == m)
    
    while (!converged && step < maxStepsDown) {
        step <- step + 1
        
        ## backup
        lambda0 <- lambda
        thr0 <- thr
        R10 <- R1
        
        stopifnot(length(R1) > 0L)
        stopifnot(length(R1) < m)
        
        mat1 <- mat[-R1, ]
        
        ## re-calibration of lambda, *holding sk fixed*
        kMax <- min(kMax, nrow(mat1))
        C <- setdiff(1:m, R1)
        pivStat <-  pivStatFUN(mat, kMax, C)
        lambda <- stats::quantile(pivStat, alpha, type = 1)
        thr <- tk(lambda)
        
        thr1 <- thr[1]   ## (1-)FWER threshold
        R1 <- which(p.values <= thr1)
        
        ## convergence reached?
        noNewRejection <- all(R1 %in% R10)
        ## In rare situations R1 is strictly included in R10,
        ## we need to declare convergence in such cases as well
        ## (and keep the largest rejection set!)
        if (noNewRejection) {
            if (!identical(R1, R10)) {
                print("strict inclusion!")
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
    # shouldn't be a warning
    # if (step == maxStepsDown && maxStepsDown > 0) {
    #     warning("Maximal number of steps down reached without reaching convergence")
    # }
    if (!is.null(p.values)) {
        ## upper bound on the number of false positives among first 'natural' rejections
        o <- order(p.values)
        Vbar <- curveMaxFP(p.values[o], thr)
    } else {
        Vbar <- NULL
    }
    stepsDown <- list(
        thr = thrMat,
        pivStat = pivMat,
        lambda = lambdas)
    
    res <- list(thr = thr,
                pivStat = pivStat,
                lambda = lambda,
                Vbar = Vbar,
                stepsDown = stepsDown)
    return(res);
}


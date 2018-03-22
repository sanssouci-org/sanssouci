#' Upper bound for the number of false discoveries in a selection
#' 
#' @param stat A vector of test statistics of the selected items
#' @param thr A vector of non-increasing \eqn{K} JER-controlling thresholds
#' @return An upper bound on the number of false discoveries in the selection
#' @export
#' @examples
#' 
#' m <- 200
#' B <- 1e3
#' rho <- 0.3
#' testStat <- gaussianTestStatistics(m, B, pi0 = 0.5, SNR = 2, dep = "equi", param = rho)
#' mat <- testStat$X0
#' stat <- testStat$x
#' cal <- jointFWERControl(mat = mat, refFamily = "Simes", alpha = 0.2)
#' cal$lambda  ## > alpha if rho > 0
#' 
#' thr <- sort(cal$thr, decreasing = TRUE)
#' 
#' M0 <- FP(stat, thr) ## upper bound on m0...
#' M0/m
#' 
#' sstat <- sort(stat, decreasing=TRUE)
#' FP(head(sstat), thr)
#' FP(c(head(sstat), tail(sstat)), thr)
#' 
FP <- function(stat, thr) {
    stopifnot(identical(sort(thr, decreasing=TRUE), thr))
    nS <- length(stat)
    K <- length(thr)
    
    size <- min(nS, K)
    seqK <- 1:size
    thr <- thr[seqK]  ## k-FWER control for k>nS is useless (will yield bound > nS)
    
    card <- sapply(thr, FUN = function(thr) {
        sum(stat < thr)
    })
    min(nS, card + seqK - 1)
}


#' Upper bound for the number of true discoveries in a selection
#' 
#' @describeIn FP
#' @return An upper bound on the number of true discoveries in the selection
#' @export
TP <- function(stat, thr) {
    length(stat) - TP(stat, thr)
}

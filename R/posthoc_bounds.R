#' Upper bound for the number of false discoveries in a selection
#' 
#' @param stat A vector of test statistics of the selected items
#' @param thr A vector of non-increasing \eqn{K} JER-controlling thresholds
#' @return An upper bound on the number of false discoveries in the selection
#' @export
#' @examples
#' 
#' m <- 123
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100, 
#'                        pi0 = 0.8, SNR = 3, prob = 0.5)
#' X <- sim$X
#' cal <- calibrateJER(X, B = 1e3, alpha = 0.2, refFamily="Simes", )
#' thr <- sort(cal$thr, decreasing = TRUE)
#' stat <- cal$stat
#' 
#' M0 <- maxFP(stat, thr) ## upper bound on m0...
#' M0/m
#' 
#' sstat <- sort(stat, decreasing=TRUE)
#' maxFP(head(sstat), thr)
#' maxFP(c(head(sstat), tail(sstat)), thr)
#' 
maxFP <- function(stat, thr) {
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


#' Lower bound for the number of true discoveries in a selection
#' 
#' @inheritParams maxFP
#' @return A Lower bound on the number of true discoveries in the selection
#' @export
minTP <- function(stat, thr) {
    length(stat) - maxFP(stat, thr)
}

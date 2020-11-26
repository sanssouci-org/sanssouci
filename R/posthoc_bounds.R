#' Upper bound for the number of false discoveries in a selection
#' 
#' @param p.values A vector of p-values for the selected items
#' @param thr A vector of non-decreasing k-FWER-controlling thresholds
#' @return A post hoc upper bound on the number of false discoveries in the selection
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @examples
#' 
#' m <- 123
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100, 
#'                        pi0 = 0.8, SNR = 2.5, prob = 0.5)
#' X <- sim$X
#' cal <- calibrateJER(X, sim$categ, B = 1e3, alpha = 0.2, refFamily="Simes", )
#' thr <- sort(cal$thr)
#' pval <- sort(cal$p.values)
#' 
#' M0 <- maxFP(pval, cal$thr) ## upper bound on m0...
#' M0/m
#' 
#' maxFP(head(pval), thr)
#' maxFP(tail(pval), thr)
#' maxFP(c(head(pval), tail(pval)), thr)
#' 
maxFP <- function(p.values, thr) {
    stopifnot(identical(sort(thr), thr))
    nS <- length(p.values)
    K <- length(thr)
    
    
    size <- min(nS, K)
    if (size == 0) {
        return(0)
    }
    seqK <- seq(from = 1, to = size, by = 1)
    thr <- thr[seqK]  ## k-FWER control for k>nS is useless (will yield bound > nS)
    
    card <- sapply(thr, FUN = function(thr) {
        sum(p.values > thr)
    })
    min(nS, card + seqK - 1)
}


#' Lower bound for the number of true discoveries in a selection
#' 
#' @inheritParams maxFP
#' @return A Lower bound on the number of true discoveries in the selection
#' @export
minTP <- function(p.values, thr) {
    length(p.values) - maxFP(p.values, thr)
}

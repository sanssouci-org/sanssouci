#' (Slow) Upper bound for the number of false discoveries in a selection
#' 
#' @param p.values A vector of p-values for the selected items
#' @param thr A vector of non-decreasing k-FWER-controlling thresholds
#' @return A post hoc upper bound on the number of false discoveries in the selection
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @seealso get_max_FP
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


#' (Slow) Lower bound for the number of true discoveries in a selection
#' 
#' @inheritParams maxFP
#' @return A Lower bound on the number of true discoveries in the selection
#' @seealso get_min_TP
#' @export
minTP <- function(p.values, thr) {
    length(p.values) - maxFP(p.values, thr)
}


#' Post hoc bounds
#' 
#' get_max_FP: upper bound on the number of false discoveries
#' 
#' @rdname posthoc-bounds
#' @param p.values A vector of p-values
#' @param thr A vector of non-decreasing JER-controlling thresholds
#' @return A function \code{f} taking in argument a subset \code{S} of indices,
#'   such that \code{f(S)} returns the value of the corresponding bound.
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @examples
#'
#' m <- 543
#' pi0 <- 0.8
#' sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
#'                        pi0 = pi0, SNR = 3, prob = 0.5)
#' X <- sim$X
#' categ <- sim$categ
#' alpha <- 0.1
#' cal <- calibrate(X, categ, B = 1e2, alpha = alpha, family = "Simes")
#' thr <- cal$thr
#' pval <- rowWelchTests(X, categ)$p.value
#' 
#' S <- which(pval < 0.02)
#' 
#' get_max_FP(pval, thr)(S)
#' get_max_FDP(pval, thr)(S)
#' get_min_TP(pval, thr)(S)
#' get_min_TDP(pval, thr)(S)
#' 
#' plot(get_max_FDP(pval, thr)(1:m, envelope=TRUE))
#' 
get_max_FP <- function(p.values, thr) {
    stopifnot(length(p.values) >= length(thr)) 
    max_FP <- function(S, envelope = FALSE) {
        p.values <- sort(p.values[S])
        res <- sansSouci::curveMaxFP(p.values, thr, flavor = "BNR2016")
        if (!envelope) {
            res <- res[length(res)]
        }
        res
    }
    return(max_FP)
}

#' 
#' @description get_min_TP: lower bound on the number of true discoveries
#' 
#' @rdname posthoc-bounds
#' @export
get_min_TP <- function(p.values, thr) {
    max_FP <- get_max_FP(p.values, thr)
    min_TP <- function(S, envelope = FALSE) {
        res <- seq(along = S) - max_FP(S = S, envelope = TRUE)
        if (!envelope) {
            res <- res[length(res)]
        }
        res
    }
    return(min_TP)
}

#' @rdname posthoc-bounds
#' @description get_max_FDP: upper bound on the false discovery proportion
#' @export
get_max_FDP <- function(p.values, thr) {
    max_FP <- get_max_FP(p.values, thr)
    max_FDP <- function(S, envelope = FALSE) {
        res <- max_FP(S = S, envelope = TRUE)/seq(along = S)
        if (!envelope) {
            res <- res[length(res)]
        }
        res <- rev(cummin(rev(res)))
        res
    }
    return(max_FDP)
}

#' @description get_min_TDP: lower bound on the true discovery proportion
#' @rdname posthoc-bounds
#' @export
get_min_TDP <- function(p.values, thr) {
    max_FDP <- get_max_FDP(p.values, thr)
    min_TDP <- function(S, envelope = FALSE) {
        1 - max_FDP(S = S, envelope = envelope)
    }
    return(min_TDP)
}

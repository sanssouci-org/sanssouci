posthoc_bound <- function(p.values, S = seq_along(p.values), thr = NULL, lab = NULL, 
                          what = c("TP", "FDP"), all = FALSE) {
    if (is.null(thr)) {
        stop("Argument 'thr' must be non NULL")
    }
    if (is.null(lab)) {
        stop("Argument 'lab' must be non NULL")
    }
    s <- length(S)
    idxs <- seq_len(s)
    max_FP <- rep(NA_integer_, s)
    pS <- p.values[S]
    o <- order(pS)
    sorted_p <- pS[o]
    
    if (length(thr) == length(p.values) && all(thr %in% c(0,1))) {
        # assume 'thr' is in fact the "truth" <=> Oracle thresholds
        # then it suffices to count the number of '0' in 'thr', cumulatively
        max_FP <- cumsum(thr[o] == 0) 
    } else {
        max_FP <- curveMaxFP(sorted_p, thr)
    }
    bounds <- formatBounds(max_FP, idxs = idxs, lab = lab, what = what, all = all)
    bounds
}


formatBounds <- function(max_FP, idxs = seq_len(max_FP), lab = NULL, 
                         what = c("TP", "FDP"), all = FALSE) {
    stopifnot(length(max_FP) == length(idxs))
    max_FDP <- max_FP/idxs
    what0 <- c("FP", "TP", "FDP", "TDP")
    if (!all(what %in% what0)) {
        stop("Error in argument 'what': only the following statistics are supported: ", paste(what0, collapse = ", "))
    }
    if (length(max_FP) == 0){ # output should be empty
        mat <- matrix(NA, nrow = 0, ncol = 4)
        colnames(mat) <- c("x", "label", "stat", "bound")
        bounds <- as.data.frame(mat)
    } else {
        annot <- data.frame(x = idxs, 
                            label = lab,
                            row.names = NULL)
        boundsList <- list(
            FP = cbind(annot, stat = "FP", bound = max_FP),
            TP = cbind(annot, stat = "TP", bound = idxs - max_FP),
            FDP = cbind(annot, stat = "FDP", bound = max_FDP),
            TDP = cbind(annot, stat = "TDP", bound = 1 - max_FDP))
        boundsList <- boundsList[what]
        if (!all) {
            boundsList <- lapply(boundsList, FUN = function(x) x[length(idxs), ])
        }
        bounds <- Reduce(rbind, boundsList)
    }
    bounds
}

#' Plot confidence bound
#' 
#' @param conf_bound A data.frame or a list of data.frames as output by [predict]
#' @param xmax Right limit of the plot
#' @param cols A vector of colors of the same length as `conf_bound`
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @examples
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' res <- fit(obj, B = 100, alpha = 0.1)
#' cb <- predict(res, all = TRUE)
#' plotConfCurve(cb, xmax = 200)  ## equivalent to 'plot(res, xmax = 200)'
#' 
#' # plot two confidence curves
#' res_beta <- fit(res, B = 100, alpha = 0.1, family = "Beta", K = 20)
#' cb_beta <- predict(res_beta, all = TRUE)
#' 
#' bounds <- list("Simes"= cb, 
#'                 "Beta" = cb_beta)
#' plotConfCurve(bounds, xmax = 200)
#' 

plotConfCurve <- function(conf_bound, xmax, cols = NULL) {
    nb <- 1
    if (class(conf_bound) == "data.frame") {    # (assume) a single conf. bound
        ## do nothing!
    } else if (class(conf_bound) == "list") {          # (assume) a list of conf. bounds
        nb <- length(conf_bound)
        nms <- names(conf_bound)
        if (!is.null(nms)) {
            for (kk in seq_along(conf_bound)) {
                conf_bound[[kk]]$Template <- nms[kk]
            }
        }
        if (!is.null(cols)) {
            stopifnot(length(conf_bound) == length(cols))
        } else {
            cols <- scales::hue_pal()(length(conf_bound))
        }
        conf_bound <- Reduce(rbind, conf_bound)
        x <- NULL; rm(x); ## To please R CMD check
    } 
    if (!missing(xmax)) {
        conf_bound <- subset(conf_bound, x <= xmax) 
    }    
    
    p <- ggplot2::ggplot(conf_bound, 
                         ggplot2::aes_string(x = "x", y = "bound"))
    if (nb > 1) {
        p <- p + ggplot2::aes_string(color = "Template", linetype = "Template")
    }
    p + ggplot2::geom_line() +
        ggplot2::facet_wrap(~ stat, scales = "free_y") + 
        ggplot2::labs(x = "Number of top features selected", 
                      y = "Post hoc confidence bounds") +
        ggplot2::theme_bw() + 
        ggplot2::theme(strip.background = NULL) +
        ggplot2::scale_color_manual(values = cols) 
}    

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
#'                        pi0 = 0.8, SNR = 3, prob = 0.5)
#' X <- sim$X
#' groups <- sim$categ
#' p <- rowWelchTests(X, groups)$p.value
#' 
#' null_groups <- replicate(100, sample(groups))
#' p0 <- rowWelchTests(X, null_groups)$p.value
#' calib <- calibrate(p0, m, alpha = 0.1)
#' thr <- calib$thr
#' 
#' M0 <- maxFP(p, thr)
#' M0/m
#' 
#' sorted_p <- sort(p)
#' maxFP(head(sorted_p, 20), thr) # some signal
#' maxFP(tail(sorted_p), thr)     # no signal
#' maxFP(c(head(sorted_p), tail(sorted_p)), thr)

maxFP <- function(p.values, thr) {
    s <- length(p.values)
    if (s == 0) {
        return(0)
    }
    all_maxFP <- curveMaxFP(p.values, thr)
    all_maxFP[s]
}


maxFP_slow <- function(p.values, thr) {
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

#' Lower bound for the true discovery proportion in a selection
#' 
#' @inheritParams maxFP
#' @return Lower bound on the proportion of true discoveries in the selection
#' @export
minTDP <- function(p.values, thr) {
    1 -  maxFP(p.values, thr)/length(p.values)
}

#' Upper bound for the false discovery proportion in a selection
#' 
#' @inheritParams maxFP
#' @return Upper bound on the proportion of false discoveries in the selection
#' @export
maxFDP <- function(p.values, thr) {
    maxFP(p.values, thr)/length(p.values)
}

#' Upper bound for the number of false discoveries among most significant items
#'
#' @param p.values A vector containing s p-values
#' @param thr A vector of \eqn{K} JER-controlling thresholds

#' @return A vector of size \eqn{s} giving an joint upper confidence bound on
#'   the number of false discoveries among the \eqn{k} most significant items
#'   for all \eqn{i \in \{1\ldots s\}}.

#' @author Gilles Blanchard, Nicolas Enjalbert-Courrech, Pierre Neuvial and
#'   Etienne Roquain
#' @details The time and space complexity of this function is O(s), which is
#'   optimal since s is the length of the returned vector.

curveMaxFP <- function(p.values, thr) {
    s <- length(p.values)
    if (s == 0) {
        return(numeric(0L))
    }
    p.values <- sort(p.values)
    thr <- sort(thr)
    
    kMax <- length(thr)
    if (s < kMax){  # truncate thr to first 's' values
        seqK <- seq(from = 1, to = s, by = 1)
        thr <- thr[seqK]
    } else { # complete 'thr' to length 's' with its last value
        thr <- c(thr, rep(thr[kMax], s - kMax))
    }
    ## sanity checks
    stopifnot(length(thr) == s)
    rm(kMax)
    
    K <- rep(s, s) ## K[i] = number of k/ T[i] <= s[k]
    Z <- rep(s, s) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
    ## 'K' and 'Z' are initialized to their largest possible value (both 's')
    kk <- 1
    ii <- 1
    while ((kk <= s) && (ii <= s)) {
        if (thr[kk] >= p.values[ii]) {
            K[ii] <- kk-1
            ii <- ii+1
        } else {
            Z[kk] <- ii-1
            kk <- kk+1
        }
    }
    Vbar <- numeric(s)
    ww <- which(K > 0)
    A <- Z - (1:s) + 1
    cA <- cummax(A)[K[ww]]  # cA[i] = max_{k<K[i]} A[k]
    Vbar[ww] <- pmin(ww - cA, K[ww])
    Vbar
}

curveMaxFP_ECN <- function(p.values, thr) {
    m <- length(p.values)
    kMax <- length(thr)
    if (m < kMax){
        seqK <- seq(from = 1, to = m, by = 1)
        thr <- thr[seqK]
        kMax <- length(thr) 
    }
    if (kMax < m) {
        thr <- c(thr, rep(thr[kMax], m-kMax))
        kMax <- length(thr)
        stopifnot(kMax==m)
    }
    ## sanity checks
    ##stopifnot(length(stat)==m)
    stopifnot(identical(sort(thr), thr))
    stopifnot(identical(sort(p.values), p.values))
    
    K <- rep(kMax, m) ## K[i] = number of k/ T[i] <= s[k] = BB in 'Mein2006'
    Z <- rep(m, kMax) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
    ## 'K' and 'Z' are initialized to their largest possible value, 
    ## ie 'm' and 'kMax', respectively
    kk <- 1
    ii <- 1
    while ((kk <= kMax) && (ii <= m)) {
        if (thr[kk] >= p.values[ii]) {
            K[ii] <- kk-1
            ii <- ii+1
        } else {
            Z[kk] <- ii-1
            kk <- kk+1
        }
    }
    Vbar <- numeric(m)
    ww <- which(K>0)
    A <- Z - (1:kMax)+1
    cA <- cummax(A)[K[ww]]  # cA[i] = max_{k<K[i]} A[k]
    Vbar[ww] <- pmin(ww-cA, K[ww])
    # Vbar[ww] <- ww-cA
    # www <- which(Vbar > K & Vbar < kMax)
    # www <- which(Vbar > K)
    # Vbar[www] <- K[www]
    Vbar
}

#' Upper bound for the number of false discoveries among most significant items (older version)
#'
#' @param p.values A vector containing all \eqn{m} p-values, sorted non-decreasingly
#' @param thr A vector of \eqn{K} JER-controlling thresholds, sorted non-decreasingly
#' @param flavor The algorithm to compute the bound 'BNR2014' and
#' 'BNR2016' give identical results. Both should be slightly better
#' than 'Mein2006' (example situation?). Flavor 'BNR2016' has a
#' linear time complexity, hence it is much faster than 'Mein2006'
#' and much much faster than 'BNR2014'.
#' @return A vector of size \eqn{m} giving an joint upper confidence bound on
#'   the number of false discoveries among the \eqn{k} most significant items
#'   for all \eqn{k \in \{1\ldots m\}}.
#' @author Gilles Blanchard, Nicolas Enjalbert-Courrech, Pierre Neuvial and Etienne Roquain
#' @details These older implementations of 'curveMaxFP' are here for the purpose of testing that the current one yields consistent results.

curveMaxFP_old <- function(p.values, thr, 
                           flavor = c("BNR2016", "Mein2006", "BNR2014")) {
    flavor <- match.arg(flavor)
    m <- length(p.values)
    kMax <- length(thr)
    if (kMax < m && flavor %in% c("Mein2006", "BNR2016")) {
        thr <- c(thr, rep(thr[kMax], m-kMax))
        kMax <- length(thr)
        stopifnot(kMax==m)
    }
    if (flavor=="Mein2006") {
        ## (loose) upper bound on number of FALSE discoveries among first rejections
        R <- 1:m
        BB <- sapply(p.values[R], function(x) sum(x>thr))     ## Eqn (7) in Meinshausen 
        ## corresponds to 'K' in 'BNR2016'
        
        ## lower bound on number of TRUE discoveries among first rejections
        Sbar <- pmax(0, cummax(R-BB))
        
        ## (tighter) upper bound on number of FALSE discoveries among first rejections
        Vbar <- R-Sbar[R]
    } else if (flavor=="BNR2014") {    ## Etienne's version
        bound <- function(kk, ii) {
            (kk-1) + sum(p.values[1:ii] > thr[kk])
        }
        Vbar <- sapply(1:m, function(ii) {
            cand <- sapply(1:kMax, bound, ii)
            min(cand)
        })
    } else if (flavor=="BNR2016") {    ## Pierre's version
        ## sanity checks
        ##stopifnot(length(stat)==m)
        stopifnot(identical(sort(thr), thr))
        stopifnot(identical(sort(p.values), p.values))
        
        K <- rep(kMax, m) ## K[i] = number of k/ T[i] <= s[k] = BB in 'Mein2006'
        Z <- rep(m, kMax) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
        ## 'K' and 'Z' are initialized to their largest possible value, 
        ## ie 'm' and 'kMax', respectively
        kk <- 1
        ii <- 1
        while ((kk <= kMax) && (ii <= m)) {
            if (thr[kk] >= p.values[ii]) {
                K[ii] <- kk-1
                ii <- ii+1
            } else {
                Z[kk] <- ii-1
                kk <- kk+1
            }
        }
        Vbar <- numeric(m)
        ww <- which(K>0)
        A <- Z - (1:kMax)+1
        cA <- cummax(A)[K[ww]]  # cA[i] = max_{k<K[i]} A[k]
        Vbar[ww] <- pmin(ww-cA, K[ww])
        # Vbar[ww] <- ww-cA
        # www <- which(Vbar > K & Vbar < kMax)
        # www <- which(Vbar > K)
        # Vbar[www] <- K[www]
    }
    Vbar
}
############################################################################
## HISTORY of 'curveMaxFP'
##
## 2016-07-04
## o Implemented the case 'kMax<m' for all flavors. Added
##   corresponding 'testthat'scripts.
##
## 2016-06-01
## o Added flavor 'BNR2016', a much faster version of flavor
## 'BNR2014'.
##
## 2014-05-22
## o Created.
############################################################################


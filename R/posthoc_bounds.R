#' Post hoc confidence bounds on the true/false positives
#' 
#' @param object An object. See individual methods for specifics
#' @param S A subset of indices
#' @param ... Further arguments to be passed to other methods
#' @export
bound <- function(object, S, ...) UseMethod("bound")

#' @rdname bound
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#' 
#' - FP: Upper bound on the number of false positives in the 'x' most significant items
#' - TP: Lower bound on the number of true positives in the 'x' most significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most significant items.
#' 
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S` be returned? If `FALSE` (the default), only the bound for `S` is returned.
#' 
#' @return If `all` is `FALSE` (the default), only the value of the bound is returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
#' - x: Number of most significant items selected
#' - label: Label for the procedure, typically of the form 'family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`.
#' 
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' res <- fit(obj, B = 100, alpha = 0.1)
#' 
#' # post hoc bounds for all hypotheses
#' bound(res)
#'
#' # confidence curve
#' cb <- bound(res, all = TRUE)
#' head(cb)
#' plot(res)
#' 
#' # confidence curve for a subset
#' S <- which(foldChanges(res) > 0.3)
#' plot(res, S = S)
#' 
#' # plot two confidence curves
#' res_beta <- fit(obj, B = 100, alpha = 0.1, family = "Beta", K = 20)
#' cb_beta <- bound(res_beta, all = TRUE)
#' 
#' bounds <- list("Simes"= cb, 
#'                    "Beta" = cb_beta)
#' plotConfCurve(bounds, xmax = 200)
#' 
#' @export
bound.SansSouci <- function(object, S = seq_len(nHyp(object)), 
                            what = c("TP", "FDP"), all = FALSE, ...) {
    p.values <- pValues(object)
    thr <- thresholds(object)
    lab <- label(object)
    if (max(S) > nHyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    bounds <- bound(object = p.values, S = S, thr = thr, lab = lab, what = what, all = all)
    if (!all) {
        bounds <- bounds[, "bound"]
        if (length(bounds) > 1) {
            names(bounds) <- what
        }
    }
    return(bounds)
}

bound.numeric <- function(object, S = seq_along(object), thr = NULL, lab = NULL, 
                     what = c("TP", "FDP"), all = FALSE) {
    if (is.null(thr)) {
        stop("Argument 'thr' must be non NULL")
    }
    p.values <- object; rm(object);
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
        max_FP <- curveMaxFP(sorted_p, thr) ## Would be faster to do 'thr[length(S)]' here. Is it correct?
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
    Reduce(rbind, boundsList)
}

#' Plot confidence bound
#' 
#' @param conf_bound A data.frame or a list of data.frames as output by 
#'   \code{\link{bound}}
#' @param xmax Right limit of the plot
#' @param cols A vector of colors of the same length as `conf_bound`
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#'
#' @export
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

#' (Slow) Upper bound for the number of false discoveries in a selection
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


#' (Slow) Lower bound for the number of true discoveries in a selection
#' 
#' @inheritParams maxFP
#' @return A Lower bound on the number of true discoveries in the selection
#' @export
minTP <- function(p.values, thr) {
    length(p.values) - maxFP(p.values, thr)
}


#' Upper bound for the number of false discoveries among most significant items
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
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain

curveMaxFP <- function(p.values, thr, flavor=c("BNR2016", "Mein2006", "BNR2014")) {
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


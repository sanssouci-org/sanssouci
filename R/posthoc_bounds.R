#' Post hoc confidence bounds on the true/false positives
#' 
#' @inheritParams nHyp
#' @param S A subset of indices
#' @export
bound <- function(object, S, ...) UseMethod("bound")

#' Post hoc confidence bounds on the true/false positives
#' 
#' @inheritParams bound
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#' 
#' \describe{
#' \item{FP}{Upper bound on the number of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the number of true positives in the 'x' most significant items}
#' \item{FDP}{Upper bound on the proportion of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the proportion of true positives in the 'x' most significant items}}
#' Defaults to \code{c("TP", "FDP")}.
#' 
#' @return If \code{all} is \code{TRUE}, a \code{data.frame} with \eqn{m} rows and 4 columns: \describe{
#' \item{x}{Number of most significant items selected}
#' \item{label}{Label for the procedure, typically of the form '<refFamily>(<param>)'}
#' \item{bound}{Value of the post hoc bound}
#' \item{stat}{Type of post hoc bound, as specified by argument \code{bound}}. If \code{all} is \code{FALSE}, only the value of the bound is returned. 
#' }
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' obj <- SansSouci(Y = sim$X, groups = sim$categ)
#' res <- fit(obj, B = 100, alpha = 0.1)
#' # calculate, print, and plot confidence bound
#' ce <- bound(res, all = TRUE)
#' head(ce)
#' 
#' ce <- bound(res, S=1:10, all = TRUE)
#' head(ce)
#' 
#' plot(res, S=which(sim$H==1))
#' 
#' @export
bound.SansSouci <- function(object, S = 1:nHyp(object), 
                            what = c("TP", "FDP"), all = FALSE) {
    p.values <- pValues(object)
    thr <- thresholds(object)
    lab <- label(object)
    if (max(S) > nHyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    bounds <- bound(object = p.values, S = S, thr = thr, lab = lab, what = what, all = all)
    if (!all) {
        bounds <- bounds[, "bound"]
        names(bounds) <- what
    }
    return(bounds)
}

bound.numeric <- function(object, S, thr, lab, 
                     what = c("TP", "FDP"), all = FALSE) {
    p.values <- object; rm(object);
    s <- length(S)
    o <- order(p.values)
    idxs <- seq_len(s)
    maxFP <- rep(NA_integer_, s)
    if (length(thr) == length(p.values) && all(thr %in% c(0,1))) {
        # assume 'thr' is in fact the "truth" <=> Oracle thresholds
        # then it suffices to count the number of '0' in 'thr', cumulatively
        maxFP <- cumsum(thr[o[S]]==0) 
    } else {
        sorted_p <- p.values[o[S]]
        maxFP <- curveMaxFP(sorted_p, thr) ## Would be faster to do 'thr[length(S)]' here. Is it correct?
        rm(sorted_p, o, p.values)
    }
    maxFDP <- maxFP/idxs
    what0 <- c("FP", "TP", "FDP", "TDP")
    if (!all(what %in% what0)) {
        stop("Error in argument 'what': only the following statistics are supported: ", paste(what0, collapse = ", "))
    }
    annot <- data.frame(x = idxs, 
                        label = lab,
                        row.names = NULL)
    boundsList <- list(
        FP = cbind(annot, stat = "FP", bound = maxFP),
        TP = cbind(annot, stat = "TP", bound = idxs - maxFP),
        FDP = cbind(annot, stat = "FDP", bound = maxFDP),
        TDP = cbind(annot, stat = "TDP", bound = 1 - maxFDP))
    boundsList <- boundsList[what]
    if (!all) {
        boundsList <- lapply(boundsList, FUN = function(x) x[s, ])
    }
    bounds <- Reduce(rbind, boundsList)
    return(bounds)
}

#' Plot confidence bound
#' 
#' @param conf_bound A data.frame or a list of data.frames as output by 
#'   \code{\link{bound}}
#'
#' @param xmax Right limit of the plot
#' @param cols A vector of colors of the same length as `conf_bound`
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#'
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.3, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' dat <- sim$X
#' categ <- sim$categ
#' rwt <- rowWelchTests(dat, categ, alternative = "greater")
#' 
#' # calculate and plot confidence bound
#' alpha <- 0.1
#' ce <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = alpha)
#' plotConfCurve(ce, xmax = 200) 
#' 
#' # calculate and plot several confidence bounds
#' B <- 100
#' cal <- calibrateJER(X = dat, categ = categ, B = B, alpha = alpha, refFamily = "Simes")
#' cal_beta <- calibrateJER(X = dat, categ = categ, B = B, alpha = alpha, refFamily = "Beta", K = 20)
#' cec <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = cal$lambda)

#' all_bounds <- list("Simes" = ce, 
#'                 "Simes + calibration"= cal$conf_bound, 
#'                 "Beta + calibration" = cal_beta$conf_bound)
#' plotConfCurve(all_bounds, xmax = 200)
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

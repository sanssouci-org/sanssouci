#'
#' Confidence envelope on the true/false positives among most significant items
#'
#' @param p.values A vector containing all \eqn{m} p-values, sorted
#'   increasingly
#'   
#' @param refFamily A character value, the reference family to be used. Should
#'   be either "Simes" (or equivalenlty, "Linear"), "Beta", or "Oracle".
#'   
#' @param param A numeric value or vector of parameters for the reference family. 
#' 
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}. Defaults to \eqn{m}.
#' 
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
#' @details \code{param} should be a numeric value unless \code{refFamily ==
#'   "Oracle"}. In the latter case, \code{param} should be a boolean vector of
#'   length \eqn{m} indicating whether each null hypothesis is true or false.
#'
#' @return A \code{data.frame} with \eqn{m} rows and 5 columns: \describe{
#' \item{x}{Number of most significant items selected}
#' \item{family}{Matches input argument \code{refFamily}}
#' \item{param}{Matches argument \code{param}}
#' \item{procedure}{Label for the procedure, typically of the form '<refFamily>(<param>)'}
#' \item{bound}{Value of the post hoc bound}
#' \item{stat}{Type of post hoc bound, as specified by argument \code{bound}}
#' }
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' rwt <- rowWelchTests(sim$X, sim$categ, alternative = "greater")
#' 
#' # calculate, print, and plot confidence envelope
#' ce <- confidenceEnvelope(rwt$p.value, refFamily = "Simes", param = 0.1)
#' head(ce)
#' plotConfidenceEnvelope(ce, xmax = 200) 
#' 
confidenceEnvelope <- function(p.values, refFamily, param, K = length(p.values), what = c("TP", "FDP")) {
    m <- length(p.values)
    idxs <- 1:m
    o <- order(p.values)
    rk <- rank(p.values)
    fam0 <- c("Simes", "Beta", "Oracle")
    if (!(refFamily %in% fam0)) {
        stop("Unknown family: ", refFamily, "\n",
             "Only the following reference families are currently supported: ", 
             paste(fam0, collapse = ", "))
    }
    what0 <- c("FP", "TP", "FDP", "TDP")
    if (!all(what %in% what0)) {
        stop("Error in argument 'what': only the following statistics are supported: ", what0)
    }
    max_FP <- rep(NA_real_, m)
    if (refFamily %in% c("Simes", "Linear")) {
        thr <- SimesThresholdFamily(m, kMax = K)(param)
        max_FP <- curveMaxFP(p.values[o], thr)
    } else if (refFamily == "Beta") {
        thr <- BetaThresholdFamily(m, kMax = K)(param)
        max_FP <- curveMaxFP(p.values[o], thr)
    } else if (refFamily == "Oracle") {
        stopifnot(length(param) == m)
        max_FP <- cumsum(o %in% which(param))
    }
    proc <- sprintf("%s(%s)", refFamily, param)
    if (refFamily == "Oracle") {
        proc <- "Oracle"
        param <- NA_real_
    }
    max_FDP <- max_FP/idxs
    annot <- data.frame(x = idxs, 
                        rank = rk,
                        family = refFamily,
                        param = param,
                        procedure = proc,
                        row.names = NULL)
    boundsList <- list(
        FP = cbind(annot, stat = "FP", bound = max_FP),
        TP = cbind(annot, stat = "TP", bound = idxs - max_FP),
        FDP = cbind(annot, stat = "FDP", bound = max_FDP),
        TDP = cbind(annot, stat = "TDP", bound = 1 - max_FDP))
    boundsList <- boundsList[what]
    Reduce(rbind, boundsList)
}

#' Plot confidence envelope
#' 
#' @param conf_env A data.frame or a list of data.frames as output by 
#'   \code{\link{confidenceEnvelope}}
#'
#' @param xmax Right limit of the plot
#' @param cols A vector of colors of the same length as `conf_env`
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
#' # calculate and plot confidence envelope
#' alpha <- 0.1
#' ce <- confidenceEnvelope(rwt$p.value, refFamily = "Simes", param = alpha)
#' plotConfidenceEnvelope(ce, xmax = 200) 
#' 
#' # calculate and plot several confidence envelopes
#' B <- 100
#' cal <- calibrateJER(X = dat, categ = categ, B = B, alpha = alpha, refFamily = "Simes")
#' cal_beta <- calibrateJER(X = dat, categ = categ, B = B, alpha = alpha, refFamily = "Beta", K = 20)
#' cec <- confidenceEnvelope(rwt$p.value, refFamily = "Simes", param = cal$lambda)

#' all_env <- list("Simes" = ce, 
#'                 "Simes + calibration"= cal$conf_env, 
#'                 "Beta + calibration" = cal_beta$conf_env)
#' plotConfidenceEnvelope(all_env, xmax = 200)
#' 
plotConfidenceEnvelope <- function(conf_env, xmax, cols = NULL) {
    nb_env <- 1
    if (class(conf_env) == "data.frame") {    # (assume) a single conf. envelope
        ## do nothing!
    } else if (class(conf_env) == "list") {          # (assume) a list of conf. envelopes
        nb_env <- length(conf_env)
        nms <- names(conf_env)
        if (!is.null(nms)) {
            for (kk in seq_along(conf_env)) {
                conf_env[[kk]]$Template <- nms[kk]
            }
        }
        if (!is.null(cols)) {
            stopifnot(length(conf_env) == length(cols))
        } else {
            cols <- scales::hue_pal()(length(conf_env))
        }
        conf_env <- Reduce(rbind, conf_env)
        x <- NULL; rm(x); ## To please R CMD check
    } 
    if (!missing(xmax)) {
        conf_env <- subset(conf_env, x <= xmax) 
    }    
    
    p <- ggplot2::ggplot(conf_env, 
                         ggplot2::aes_string(x = "x", y = "bound"))
    if (nb_env > 1) {
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

#'
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
## HISTORY:
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


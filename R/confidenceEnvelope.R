#'
#' Confidence envelope on the true/false positives among most significant items
#'
#' 
#'
#' @param stat A vector containing all \eqn{m} test statistics, sorted
#'   non-increasingly
#' @param alphas A numeric vector of target confidence level(s)
#' 
#' @param family A character vector, the name of the reference families to be
#'   used. Currently, only the "Simes" (a.k.a. "linear") and "Beta" families 
#'   are implemented.
#'   
#' @param what A character vector, the name of the statistics to be computed,
#'   among: \describe{
#' \item{FP}{Upper bound on the number of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the number of true positives in the 'x' most significant items}
#' \item{FDP}{Upper bound on the proportion of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the proportion of true positives in the 'x' most significant items}}
#'   
#'
#' @return A \code{data.frame} with \eqn{m} rows and 5 columns:\describe{
#' \item{x}{Number of top items in the selection}
#' \item{stat}{Name of the statistic, corresponding to argument \code{what}}
#' \item{bound}{Value of the post hoc bound}
#' \item{family}{Name of the reference family}
#' \item{alpha}{Target confidence level}
#' }
#' For example, a row with \code{x=10}, \code{stat="TP"}, \code{bound = 4}, \code{family="Simes"}, \code{alpha = 0.1} means that the post hoc bound derived from the Simes reference family at level 0.1 implies that at least 4 the 10 most significant items is are true positives.
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#' m <- 123
#' alpha <- 0.1
#' sim <- gaussianSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' cal <- calibrateJER(sim$X, B = 1e3, alpha = alpha, refFamily="Simes")
#' conf_env <- confidenceEnvelope(cal$stat, refFamily="Simes", alpha = cal$lambda)
#' plot(1:m, conf_env$FP, t = 's', col = "purple")
#'
#' # compare with true number of false positives
#' o <- order(cal$stat, decreasing = TRUE)
#' H0 <- which(sim$H == 0)
#' V <- cumsum(o %in% H0)
#' lines(1:m, V, col = 2)
#'
#' # compare with Simes bound (ie with no calibration)
#' conf_env_Simes <- confidenceEnvelope(cal$stat, refFamily="Simes", alpha = alpha)
#' lines(1:m, conf_env_Simes$FP, col = 1)

confidenceEnvelope <- function(stat, family, alpha, what = c("FP", "TP", "FDP", "TDP")) {
    fam0 <- c("Simes", "Beta")
    if (!all(family %in% fam0)) {
        stop("Only the following reference families are currently supported: ", fam0)
    }
    if (min(alpha) < 0 || max(alpha) > 1) {
        stop("Target level 'alpha' should be in [0,1]")
    }
    what0 <- c("FP", "TP", "FDP", "TDP")
    if (!all(what %in% what0)) {
        stop("Only the following statistics are supported: ", what0)
    }
    
    m <- length(stat)
    idxs <- 1:m
    o <- order(stat, decreasing = TRUE)
    configs <- expand.grid(family = family, alpha = alpha)
    res <- NULL
    for (kk in 1:nrow(configs)) {
        fam <- configs[kk, "family"]
        al <- configs[kk, "alpha"]
        if (fam == "Simes") {
            thr <- SimesThresholdFamily(m)(al)
        } else if (fam == "Beta") {
            thr <- BetaThresholdFamily(m)(al)
        }
        max_FP <- curveMaxFP(stat[o], thr)
        max_FDP <- max_FP/idxs
        bounds <- list(FP = max_FP,
                       TP = idxs - max_FP, 
                       FDP = max_FDP,
                       TDP = 1 - max_FDP)
        
        res_fam <- NULL
        for (ww in what) {
            res_fam <- rbind(res_fam, 
                             data.frame(x = idxs, stat = ww, bound = bounds[[ww]]))
        }
        res_fam <- cbind(res_fam, 
                         family = fam,
                         alpha = al, 
                         row.names = NULL)
        res <- rbind(res, res_fam)
    }
    res
}

#'
#' Upper bound for the number of false discoveries among most significant items
#'
#' @param stat A vector containing all \eqn{m} test statistics, sorted non-increasingly
#' @param thr A vector of \eqn{K} JER-controlling thresholds, sorted non-increasingly
#' @param flavor The algorithm to compute the bound 'BNR2014' and
#' 'BNR2016' give identical results. Both should be slightly better
#' than 'Mein2006' (example situation?). Flavor 'BNR2016' has a
#' linear time complexity, hence it is much faster than 'Mein2006'
#' and much much faster than 'BNR2014'.
#' @return A vector of size \eqn{m} giving an joint upper confidence bound on
#'   the number of false discoveries among the \eqn{k} most significant items
#'   for all \eqn{k \in \{1\ldots m\}}.
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @examples
#' m <- 123
#' alpha <- 0.2
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100, pi0 = 0.5, SNR = 3, prob = 0.5)
#' cal <- calibrateJER(sim$X, B = 1e3, alpha = alpha, refFamily="Simes")
#' stat <- cal$stat
#' o <- order(stat, decreasing = TRUE)
#' ub <- curveMaxFP(stat[o], cal$thr)
#' plot(1:m, ub, t = 's', col = "purple")
#' 
#' # compare with true number of false positives
#' H0 <- which(sim$H == 0)
#' V <- cumsum(o %in% H0)
#' lines(1:m, V, col = 2)
#' 
#' # compare with Simes bound
#' thrSimes <- SimesThresholdFamily(m)(alpha)
#' ubSimes <- curveMaxFP(stat[o], thrSimes)
#' lines(1:m, ubSimes, col = 1)

curveMaxFP <- function(stat, thr, flavor=c("BNR2016", "Mein2006", "BNR2014")) {
    m <- length(stat)
    kMax <- length(thr)

    flavor <- match.arg(flavor)
    if (flavor=="Mein2006") {
        ## (loose) upper bound on number of FALSE discoveries among first rejections
        R <- 1:m
        BB <- sapply(stat[R], function(x) sum(x<=thr))     ## Eqn (7) in Meinshausen
        stopifnot(all(BB<=kMax))  ## sanity check

        ## lower bound on number of TRUE discoveries among first rejections
        Sbar <- pmax(0, cummax(R-BB))

        ## (tighter) upper bound on number of FALSE discoveries among first rejections
        Vbar <- R-Sbar[R]
    } else if (flavor=="BNR2014") {    ## Etienne's version
        bound <- function(kk, ii) {
            (kk-1) + sum(stat[1:ii] <= thr[kk])
        }
        Vbar <- sapply(1:m, function(ii) {
                           cand <- sapply(1:kMax, bound, ii)
                           min(cand)
                       })
    } else if (flavor=="BNR2016") {    ## Pierre's version
        ## sanity checks
        ##stopifnot(length(stat)==m)
        stopifnot(identical(sort(thr, decreasing=TRUE), thr))
        stopifnot(identical(sort(stat, decreasing=TRUE), stat))
        
        K <- rep(kMax, m) ## K[i] = number of k/ T[i] <= s[k]
        Z <- rep(m, kMax) ## Z[k] = number of i/ T[i] >  s[k]
        ## 'K' and 'Z' are initialized to their largest possible value, 
        ## ie 'm' and 'kMax', respectively
        kk <- 1
        ii <- 1
        while ((kk <= kMax) && (ii <= m)) {
            if (thr[kk]<=stat[ii]) {
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
        cA <- cummax(A)[K[ww]]
        Vbar[ww] <- pmin(ww-cA, K[ww])
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


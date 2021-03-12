#' Confidence bounds on the true/false positives among most significant items
#'
#' @param p.values A vector containing the p-values for all m hypotheses, sorted
#'   increasingly
#'   
#' @param refFamily A character value, the reference family to be used. Should
#'   be either "Simes" (or equivalenlty, "Linear"), "Beta", or "Oracle".
#'   
#' @param param A numeric value or vector of parameters for the reference family. 
#' 
#' @param K For JER control over `1:K`, ie joint control of all
#'   k-FWER, k <= K. Defaults to m.
#' 
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#' 
#' - FP: Upper bound on the number of false positives in the 'x' most significant items
#' - TP: Lower bound on the number of true positives in the 'x' most significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most significant items.
#' 
#' Defaults to `c("TP", "FDP")`
#' 
#' @details `param` should be a numeric value unless `refFamily ==
#'   "Oracle"`. In the latter case, `param`` should be a boolean vector of
#'   length m indicating whether each null hypothesis is true or false.
#'
#' @return A `data.frame` with `m` rows and 5 columns:
#' - x: Number of most significant items selected
#' - family: Matches input argument `refFamily`
#' - param: Matches argument `param`
#' - procedure: Label for the procedure, typically of the form 'refFamily(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`
#' 

#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' rwt <- rowWelchTests(sim$X, sim$categ, alternative = "greater")
#' 
#' # calculate, print, and plot confidence bound
#' cb <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = 0.1)
#' head(cb)
#' plotConfCurve(cb, xmax = 200) 
#' 
confCurveFromFam <- function(p.values, refFamily, param, K = length(p.values), what = c("TP", "FDP")) {
    m <- length(p.values)
    fam0 <- c("Simes", "Beta", "Oracle")
    if (!(refFamily %in% fam0)) {
        stop("Unknown family: ", refFamily, "\n",
             "Only the following reference families are currently supported: ", 
             paste(fam0, collapse = ", "))
    }
    thr <- NULL
    if (refFamily %in% c("Simes", "Linear")) {
        thr <- SimesThresholdFamily(m, kMax = K)(param)
    } else if (refFamily == "Beta") {
        thr <- BetaThresholdFamily(m, kMax = K)(param)
    } else if (refFamily == "Oracle") {
        stopifnot(length(param) == m && all(param %in% c(0,1)))
        thr <- param
    }
    proc <- sprintf("%s(%s)", refFamily, param)
    bound(p.values, S = 1:m, thr, lab = proc, what = what, all = TRUE)
}


<<<<<<< HEAD
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
#' @export
#' 
#' @examples
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pval <-  rowWelchTests(X, categ)$p.value
#' cal <- calibrate(X, categ, B = 1e2, alpha = 0.1, family = "Simes")
#' 
#' curveMaxFP(sort(pval), cal$thr)
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
=======

>>>>>>> develop

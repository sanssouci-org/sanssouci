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
#' @export
#' @example
#' m <- 123
#' alpha <- 0.2
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100, pi0 = 0.5, SNR = 3, prob = 0.5)
#' cal <- calibrateJER(sim$X, B = 1e3, alpha = alpha, refFamily="Simes")
#' stat <- cal$stat
#' o <- order(stat, decreasing = TRUE)
#' ub <- curveMaxFP(sstat, cal$thr)
#' plot(1:m, ub, t = 's', col = "purple")
#' 
#' # compare with true number of false positives
#' H0 <- which(sim$H == 0)
#' V <- cumsum(o %in% H0)
#' lines(1:m, V, col = 2)
#' 
#' # compare with Simes bound
#' thrSimes <- SimesThresholdFamily(m)(alpha)
#' ubSimes <- curveMaxFP(sstat, thrSimes)
#' lines(1:m, ubSimes, col = 1)

curveMaxFP <- function(stat, thr, flavor=c("BNR2016", "Mein2006", "BNR2014")) {
    m <- length(stat)
    kMax <- length(thr)

    ## sanity checks
    ##stopifnot(length(stat)==m)
    stopifnot(identical(sort(thr, decreasing=TRUE), thr))
    stopifnot(identical(sort(stat, decreasing=TRUE), stat))

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
        K <- rep(kMax, m) ## K[i] = number of k/ T[i] <= s[k]
        Z <- rep(m, kMax) ## Z[k] = number of i/ T[i] >  s[k]
        ## 'K' and 'Z' are initialized to their largest possible value, ie 'm' and 'kMax', respectively
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


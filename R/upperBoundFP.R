##' Upper bound for the number of false discoveries
##'
##'
##'
##' @param stat ordered test statistics
##' @param thr JFWER thresholds
##' @param flavor The algorithm to compute the bound 'BNR2014' and
##' 'BNR2016' give identical resutls and should be slightly better
##' than 'Mein2006' (example situation?). Flavor 'BNR2016' has a
##' linear time complexity, hence it is much faster than 'Mein2006'
##' and much much faster than 'BNR2014'.
##' @return An upper bound on the number of false discoveries
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export


upperBoundFP <- function(stat, thr, flavor=c("BNR2016", "Mein2006", "BNR2014")) {
    m <- length(thr)

  ## sanity checks
  stopifnot(length(stat)==m)
  stopifnot(identical(sort(thr, decreasing=TRUE), thr))
  stopifnot(identical(sort(stat, decreasing=TRUE), stat))

  flavor <- match.arg(flavor)
  if (flavor=="Mein2006") {
    ## (loose) upper bound on number of FALSE discoveries among first rejections
    BB <- sapply(stat, function(x) sum(x<=thr))     ## Eqn (7) in Meinshausen
    R <- 1:m

    ## lower bound on number of TRUE discoveries among first rejections
    Sbar <- pmax(0, cummax(R-BB[R]))

    ## (tighter) upper bound on number of FALSE discoveries among first rejections
    Vbar <- R-Sbar[R]
  } else if (flavor=="BNR2014") {    ## Etienne's version
    bound <- function(kk, ii) {
      (kk-1) + sum(stat[1:ii] <= thr[kk])
    }
    Vbar <- sapply(1:m, function(ii) {
      cand <- sapply(1:m, bound, ii)
      min(cand)
    })
  } else if (flavor=="BNR2016") {    ## Pierre's version
    K <- rep(m, m) ## K[i] = number of k/ T[i] <= s[k]
    Z <- rep(m, m) ## Z[k] = number of i/ T[i] >  s[k]
    ## both K and Z are initialized to 'm'
    kk <- 1
    ii <- 1
    while ((kk <= m) && (ii <= m)) {
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
    A <- Z - (1:m)+1
    cA <- cummax(A)[K[ww]]
    Vbar[ww] <- pmin(ww-cA, K[ww])
  }
  Vbar
}
############################################################################
## HISTORY:
## 2016-06-01
## o Added flavor 'BNR2016', a much faster version of flavor
## 'BNR2014'.
## 2014-05-22
## o Created.
############################################################################


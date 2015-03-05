upperBoundFP <- function(### Upper bound for the number of false discoveries
    stat,
### ordered test statistics
    thr,
### JFWER thresholds
    flavor=c("Roquain2014", "Mein2006")
### 'Roquain2014' should be slightly better
    ) {
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
  } else if (flavor=="Roquain2014") {    ## Etienne's version (!! slow !!)
    bound <- function(kk, ii) {
      (kk-1) + sum(stat[1:ii] <= thr[kk])
    }
    Vbar <- sapply(1:m, function(ii) {
      cand <- sapply(1:m, bound, ii)
      min(cand)
    })
  }
  Vbar  ## A bound on the number of false discoveries
}
############################################################################
## HISTORY:
## 2014-05-22
## o Created.
############################################################################


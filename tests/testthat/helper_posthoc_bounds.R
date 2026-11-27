# Older implementations of 'curve_max_fp' kept for testing purposes
#
# Naive implementation of interpolation bound (O(m^2))
curve_max_fp_bnr2014 <- function(p_values, thr) {
  m <- length(p_values)
  k_max <- length(thr)
  bound <- function(kk, ii) {
    (kk - 1) + sum(p_values[1:ii] > thr[kk])
  }
  V_BAR <- sapply(1:m, function(ii) {
    cand <- sapply(1:k_max, bound, ii)
    min(cand)
  })
  V_BAR
}

# Confidence envelope of Meinshausen 2006
curve_max_fp_mein2006 <- function(p_values, thr) {
  m <- length(p_values)
  k_max <- length(thr)
  if (k_max < m) {
    thr <- c(thr, rep(thr[k_max], m - k_max))
    k_max <- length(thr)
  }

  ## (loose) upper bound on number of FALSE discoveries among first rejections
  R <- 1:m
  BB <- sapply(p_values[R], function(x) sum(x > thr)) ## Eqn (7) in Meinshausen
  ## corresponds to 'K' in 'BNR2016'

  ## lower bound on number of TRUE discoveries among first rejections
  S_BAR <- pmax(0, cummax(R - BB))

  ## (tighter) upper bound on number of FALSE discoveries among first rejections
  V_BAR <- R - S_BAR[R]
  V_BAR
}

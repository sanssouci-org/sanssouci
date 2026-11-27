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

# Naive confidence envelope computations on a pathway in the structured setting
# Repeatedly calls [pred_subset_structured()] on each \eqn{S_t},
# which is not optimized and time-consuming,
# this should not be used in practice and is kept for testing purposes.
pred_pathway_structured_naive <- function(
  perm,
  forest_structure,
  zetas_list,
  leaf_list,
  pruning = FALSE,
  delete_gaps = FALSE
) {
  vstars <- numeric(length(perm))

  # the naive version doesn't need a proper completion of the
  # forest structure because pred_subset_structured
  # implicitly completes, and for the same reason
  # it is allowed to prune leaves

  if (pruning) {
    pruned <- pruning(
      forest_structure,
      zetas_list,
      leaf_list,
      prune_leafs = TRUE,
      delete_gaps = delete_gaps
    )
    forest_structure <- pruned$forest_structure
    zetas_list <- pruned$zetas_list
    m <- length(unlist(leaf_list))
    if (length(perm) == m) {
      # means that the last bound to compute
      # is V^*({1, ..., m}),
      # but the pruning already computed
      # V^*({1, ..., m}) as a by-product so we
      # might as well use it:
      vstars[m] <- pruned$v_star_all
      perm <- perm[-m]
    }
  }

  S <- numeric(0)
  j <- 0
  for (t in perm) {
    j <- j + 1
    S <- c(S, t)
    vstars[j] <- pred_subset_structured(
      S,
      forest_structure,
      zetas_list,
      leaf_list
    )
  }
  return(vstars)
}

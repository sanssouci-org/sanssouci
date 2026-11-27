#' Post hoc bounds of a selection
#'
#' @param p_values a vector of p-values
#' @param S the subset of indice of the p_values vector
#' @param thr a vector a thresholds given by [calibrate]
#' @param lab label of the method used
#' @param what type of quantities. Should be in c("FP", "TP", "FDP", "TDP").
#' @param all Boolean to specify if the function return bounds for all tests
#'
#' @returns data frame with the post hoc bounds
#' @export
#'
#' @examples
#' set.seed(0xBEEF)
#' m <- 50
#' sim <- sim_gaussian_equi(
#'   m = m, rho = 0.3, n = 45,
#'   pi0 = 0.8, SNR = 3, prob = 0.5
#' )
#' Y <- sim$X
#' groups <- sim$groups
#' p <- row_welch_tests(Y, groups)$p_value
#'
#' B <- 100
#' null_groups <- replicate(B, sample(groups))
#' p0 <- row_welch_tests(Y, null_groups)$p_value
#'
#' calib0 <- calibrate0(p0, m, alpha = 0.1)
#' thr = calib0$thr
#' res <- posthoc_bound(p, thr = thr, lab = "Simes")
posthoc_bound <- function(p_values, S = seq_along(p_values), thr = NULL,
                          lab = NULL,
                          what = c("TP", "FDP"), all = FALSE) {
  if (is.null(thr)) {
    stop("Argument 'thr' must be non NULL")
  }
  if (is.null(lab)) {
    stop("Argument 'lab' must be non NULL")
  }
  if (any(is.na(p_values)) || min(p_values) < 0 || max(p_values) > 1) {
    stop("Argument 'p_values' should only contain elements between 0 and 1")
  }
  if ((length(S) > 0) && (min(S) <= 0 || max(S) > length(p_values))) {
    stop("Argument 'S' should be a subset of indices between 1 and ",
         length(p_values))
  }

  s <- length(S)
  idxs <- seq_len(s)
  ps <- p_values[S]

  o <- order(ps)
  sorted_p <- ps[o]

  max_fp <- rep(NA_integer_, s)
  if (length(thr) == length(p_values) && all(thr %in% c(0, 1))) {
    # assume 'thr' is in fact the "truth" <=> Oracle thresholds
    # then it suffices to count the number of '0' in 'thr', cumulatively
    max_fp <- cumsum(thr[S[o]] == 0)
  } else {
    max_fp <- curve_max_fp(sorted_p, thr)
  }
  bounds <- format_bounds(max_fp, idxs = idxs, lab = lab, what = what,
                          all = all)
  return(bounds)
}


#' Table of bounds
#'
#' @param max_fp a vector of the upper boud for the number of False Discovery.
#' @param idxs vector of indices of the selection.
#' @param lab label of the method used.
#' @param what type of quantities. Should be in c("FP", "TP", "FDP", "TDP").
#' @param all Boolean to specify if the function return bounds for all tests.
#'
#' @return data frame with the post hoc bounds
#' @export
#'
format_bounds <- function(max_fp, idxs = seq_len(max_fp), lab = NULL,
                          what = c("TP", "FDP"), all = FALSE) {
  stopifnot(length(max_fp) == length(idxs))
  max_fdp <- max_fp / idxs
  what0 <- c("FP", "TP", "FDP", "TDP")
  if (!all(what %in% what0)) {
    stop("Error in argument 'what':",
         " only the following statistics are supported: ",
         paste(what0, collapse = ", "))
  }
  if (length(max_fp) == 0) {
    if (all) { # output should be empty (no subsets of positive size)
      mat <- matrix(NA, nrow = 0, ncol = 4)
      colnames(mat) <- c("x", "label", "stat", "bound")
      bounds <- as.data.frame(mat)
      return(bounds)
    } else { # output should not be empty (number of FP in empty set is 0)
      idxs <- 0
      max_fp <- 0 # number of FP in empty set is 0
      max_fdp <- 0 # FDP in empty set is 0
    }
  }
  annot <- data.frame(
    x = idxs,
    label = lab,
    row.names = NULL
  )
  bounds_list <- list(
    FP = cbind(annot, stat = "FP", bound = max_fp),
    TP = cbind(annot, stat = "TP", bound = idxs - max_fp),
    FDP = cbind(annot, stat = "FDP", bound = max_fdp),
    TDP = cbind(annot, stat = "TDP", bound = 1 - max_fdp)
  )
  bounds_list <- bounds_list[what]
  if (!all) {
    bounds_list <- lapply(bounds_list, FUN = function(x) x[length(idxs), ])
  }
  bounds <- Reduce(rbind, bounds_list)
  bounds
}

#' Plot confidence bound
#'
#' @param conf_bound A data.frame or a list of data.frames as output by
#'  [predict]
#' @param xmax Right limit of the plot
#' @param cols A vector of colors of the same length as `conf_bound`
#' @param legend_title a character for the title of the legend
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc confidence bounds on false positives using reference families.
#' Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @importFrom rlang .data
#' @examples
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
#'                     prob = 0.5)
#' res <- fit(obj, B = 100, alpha = 0.1)
#' cb <- predict(res, all = TRUE)
#' plot_conf_curve(cb, xmax = 200) ## equivalent to 'plot(res, xmax = 200)'
#'
#' # plot two confidence curves
#' res_beta <- fit(res, B = 100, alpha = 0.1, family = "Beta", K = 20)
#' cb_beta <- predict(res_beta, all = TRUE)
#'
#' bounds <- list(
#'   "Simes" = cb,
#'   "Beta" = cb_beta
#' )
#' plot_conf_curve(bounds, xmax = 200)
#'
plot_conf_curve <- function(conf_bound, xmax, cols = NULL,
                            legend_title = "Contrast") {
  nb <- 1
  if (inherits(conf_bound, "data.frame")) { # (assume) a single conf. bound
    ## do nothing!
  } else if (inherits(conf_bound, "list")) { # (assume) a list of conf. bounds
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
    x <- NULL
    rm(x) ## To please R CMD check
  }
  if (!missing(xmax)) {
    conf_bound <- subset(conf_bound, x <= xmax)
  }

  p <- ggplot2::ggplot(
    conf_bound,
    ggplot2::aes(
      x = .data[["x"]],
      y = .data[["bound"]]
    )
  )
  if (nb > 1) {
    p <- p + ggplot2::aes(
      color = .data[["Template"]],
      linetype = .data[["Template"]]
    )
  }
  p + ggplot2::geom_line() +
    ggplot2::facet_wrap(~stat, scales = "free_y") +
    ggplot2::labs(
      x = "Number of top features selected",
      y = "Post hoc confidence bounds",
      color = legend_title,
      linetype = legend_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = NULL) +
    ggplot2::scale_color_manual(values = cols)
}

#' Upper bound for the number of false discoveries in a selection
#'
#' @param p_values A vector of p-values for the selected items
#' @param thr A vector of non-decreasing JER-controlling thresholds
#' @return A post hoc upper bound on the number of false discoveries in the
#' selection
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc confidence bounds on false positives using reference families.
#' Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @examples
#'
#' m <- 123
#' sim <- sim_gaussian_equi(
#'   m = m, rho = 0.2, n = 100,
#'   pi0 = 0.8, SNR = 3, prob = 0.5
#' )
#' X <- sim$X
#' groups <- sim$groups
#' p <- row_welch_tests(X, groups)$p_value
#'
#' null_groups <- replicate(100, sample(groups))
#' p0 <- row_welch_tests(X, null_groups)$p_value
#' calib <- calibrate(p0, m, alpha = 0.1)
#' thr <- calib$thr
#'
#' M0 <- max_fp(p, thr)
#' M0 / m
#'
#' sorted_p <- sort(p)
#' max_fp(head(sorted_p, 20), thr) # some signal
#' max_fp(tail(sorted_p), thr) # no signal
#' max_fp(c(head(sorted_p), tail(sorted_p)), thr)
max_fp <- function(p_values, thr) {
  s <- length(p_values)
  if (s == 0) {
    return(0)
  }
  all_max_fp <- curve_max_fp(p_values, thr)
  all_max_fp[s]
}

#' Lower bound for the number of true discoveries in a selection
#'
#' @inheritParams max_fp
#' @return A Lower bound on the number of true discoveries in the selection
#' @export
min_tp <- function(p_values, thr) {
  length(p_values) - max_fp(p_values, thr)
}

#' Lower bound for the true discovery proportion in a selection
#'
#' @inheritParams max_fp
#' @return Lower bound on the proportion of true discoveries in the selection
#' @export
min_tdp <- function(p_values, thr) {
  1 - max_fp(p_values, thr) / length(p_values)
}

#' Upper bound for the false discovery proportion in a selection
#'
#' @inheritParams max_fp
#' @return Upper bound on the proportion of false discoveries in the selection
#' @export
max_fdp <- function(p_values, thr) {
  max_fp(p_values, thr) / length(p_values)
}

#' Upper bound for the number of false discoveries among most significant items
#'
#' @param p_values A vector containing s p-values
#' @param thr A vector of \eqn{K} JER-controlling thresholds
#'
#' @return A vector of size \eqn{s} giving an joint upper confidence bound on
#'   the number of false discoveries among the \eqn{k} most significant items
#'   for all \eqn{k \in \{1\ldots s\}}.
#'
#' @export
#'
#' @author Gilles Blanchard, Nicolas Enjalbert-Courrech, Pierre Neuvial and
#'   Etienne Roquain
#' @references Enjalbert-Courrech, N. & Neuvial, P. (2022). Powerful and
#' interpretable control of false discoveries in two-group differential
#' expression studies. Bioinformatics. doi: 10.1093/bioinformatcs/btac693
#' @details The time and space complexity of this function is O(s), which is
#'   optimal since s is the length of the returned vector.
curve_max_fp <- function(p_values, thr) {
  s <- length(p_values)
  if (s == 0) {
    return(numeric(0L))
  }
  p_values <- sort(p_values)
  thr <- sort(thr)

  k_max <- length(thr)
  if (s < k_max) { # truncate thr to first 's' values
    seq_k <- seq(from = 1, to = s, by = 1)
    thr <- thr[seq_k]
  } else { # complete 'thr' to length 's' with its last value
    thr <- c(thr, rep(thr[k_max], s - k_max))
  }
  ## sanity checks
  stopifnot(length(thr) == s)
  rm(k_max)

  K <- rep(s, s) ## K[i] = number of k/ T[i] <= s[k]
  Z <- rep(s, s) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
  ## 'K' and 'Z' are initialized to their largest possible value (both 's')
  kk <- 1
  ii <- 1
  while ((kk <= s) && (ii <= s)) {
    if (thr[kk] > p_values[ii]) {
      K[ii] <- kk - 1
      ii <- ii + 1
    } else {
      Z[kk] <- ii - 1
      kk <- kk + 1
    }
  }
  V_BAR <- numeric(s)
  ww <- which(K > 0)
  A <- Z - (1:s) + 1
  CUMMAX_A <- cummax(A)[K[ww]] # cA[i] = max_{k<K[i]} A[k]
  V_BAR[ww] <- pmin(ww - CUMMAX_A, K[ww])
  V_BAR
}


############################################################################
## HISTORY of 'curveMaxFP'
##
## 2026-11
## o Rename curveMaxFP -> curve_max_fp
##
## 2021-09
## o Implemented the case |p.values| < |thr|.
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

#' Perform JER calibration from randomization p-values
#'
#' @inheritParams get_pivotal_stat
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#' @param family A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes", or "Oracle". "Linear" and "Simes" families are
#'   identical.
#' @param p A vector of 'm' p.values, used for step-down control. If not
#'   provided, single step control is performed.
#' @param max_steps_down A numeric value, the maximum number of steps down to
#'   perform. Defaults to 10 (but the algorithm generally converges in 1 or 2
#'   steps).
#' @param piv_stat0 Don't use! Should be removed soon...
#'
#' @return A list with elements
#' - thr: A numeric vector of length K, such that the estimated probability that
#' there exists an index k between 1 and K such that the k-th maximum of the
#' test statistics of is greater than `thr[k]`, is less than alpha
#' - piv_stat: A vector of `B` pivotal statitics
#' - lambda: A numeric value, the result of the calibration
#'
#' @details 'calibrate0' performs single step calibration,
#'   whereas 'calibrate' performs step-down calibration. Hence
#'   the output of 'calibrate(..., max_steps_down = 0)' and
#'   'calibrate0(...)' should be identical
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#' @export
#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 50
#' sim <- gaussianSamples(
#'   m = m, rho = 0.3, n = 45,
#'   pi0 = 0.8, SNR = 3, prob = 0.5
#' )
#' Y <- sim$X
#' groups <- sim$categ
#' p <- rowWelchTests(Y, groups)$p.value
#'
#' B <- 100
#' null_groups <- replicate(B, sample(groups))
#' p0 <- rowWelchTests(Y, null_groups)$p.value
#'
#' calib0 <- calibrate0(p0, m, alpha = 0.1) # single step
#' calib <- calibrate(p0, m, alpha = 0.1, p = p)
#' calib$lambda >= calib0$lambda
#'
#' maxFP(p, calib$thr)
#'
#' \dontrun{
#' # Gene expression
#' data(expr_ALL, package = "sanssouci.data")
#' X <- expr_ALL
#' rm(expr_ALL)
#' groups <- ifelse(colnames(X) == "BCR/ABL", 1, 0) # map to 0/1
#'
#' null_groups <- replicate(500, sample(groups))
#' perm <- rowWelchTests(X, null_groups)
#' p0 <- perm$p.value
#'
#' alpha <- 0.1
#' m <- nrow(X)
#' p <- rowWelchTests(X, groups)$p.value
#' calib_L <- calibrate(p0, m, alpha, family = "Linear")
#' calib_B <- calibrate(p0, m, alpha, family = "Beta", K = 100)
#'
#' ## post hoc bounds
#' thr <- calib_L$thr
#' minTP(p, thr) ## lower bound on true positives
#'
#'
#' ## example of user selection: rejections of BH(0.05) procedure
#' adjp <- p.adjust(p, method = "BH")
#' sel <- which(adjp < 0.05)
#' length(sel)
#'
#' minTP(p[sel], thr)
#'
#' # confidence bound on the FDP
#' FDP_bound <- sanssouci:::curveMaxFP(sort(p), thr) / seq(along = p)
#' plot(head(FDP_bound, 300),
#'   t = "s",
#'   xlab = "Number of features",
#'   ylab = "Upper bound on False Discovery Proportion"
#' )
#' }
#'
#' @export
calibrate <- function(p0, m, alpha,
                      family = c("Linear", "Beta", "Simes"),
                      K = nrow(p0),
                      p = NULL,
                      max_steps_down = 10L,
                      piv_stat0 = NULL) {
  step <- 0
  cal <- calibrate0(p0, m, alpha,
    family = family,
    K = K,
    piv_stat0 = piv_stat0
  )
  thr <- cal$thr[1] ## (1-)FWER threshold
  R1 <- integer(0)
  if (!is.null(p)) {
    R1 <- which(p <= thr)
  }

  ## force 'convergence' if nothing to gain) or nothing left to be rejected
  converged <- (length(R1) == 0L) | (length(R1) == m)

  while (!converged && step < max_steps_down) {
    step <- step + 1

    p1 <- p0[-R1, ]
    thr <- cal$thr
    lambda <- cal$lambda
    cal <- calibrate(p1, m, alpha,
      family = family,
      K = K,
      piv_stat0 = NULL
    ) ## force piv stat calc
    R1_new <- which(p < cal$thr[1])

    noNewRejection <- all(R1_new %in% R1) ## convergence reached?
    if (noNewRejection) {
      if (!identical(R1_new, R1)) { ## can this happen?? covr check
        print("Rejecting less hypotheses at next step!")
        ## not a 'TRUE' convergence: override the last step down!
        cal$thr <- thr
        cal$lambda <- lambda
      }
      converged <- TRUE ## stop the step-down process
    } else {
      R1 <- R1_new
    }
  }
  cal$steps_down <- step
  cal
}

#' @rdname calibrate
#' @export
calibrate0 <- function(p0, m, alpha,
                       family = c("Linear", "Beta", "Simes"),
                       K = nrow(p0),
                       piv_stat0 = NULL) {
  K <- force(K)
  family <- match.arg(family)
  if (family %in% c("Linear", "Simes")) {
    t_inv <- t_inv_linear
    t_ <- t_linear
  } else if (family == "Beta") {
    t_inv <- t_inv_beta
    t_ <- t_beta
  }
  pivStat <- piv_stat0
  if (is.null(pivStat)) {
    pivStat <- get_pivotal_stat(p0, m, t_inv, min(nrow(p0), K))
  }
  lambda <- stats::quantile(pivStat, alpha, type = 1)
  thr <- t_(lambda, 1:K, m)
  res <- list(
    thr = thr,
    piv_stat = pivStat,
    lambda = lambda
  )
  return(res)
}

#' Get permutation  p-values
#'
#' Get a matrix of p-values under the null hypothesis obtained by repeated permutation of class labels (two-sample tests)
#'
#' @param X A matrix of `m` variables (hypotheses) by `n` observations
#' @param categ An numeric vector of `n` values in `0,1`
#'   specifying the column indices of the first and second samples.
#' @param B An integer value, the number of permutations to be performed
#' @param rowTestFUN a vectorized testing function (same I/O as [rowWelchTests])
#' @param alternative A character string specifying the alternative hypothesis, to be passed to `rowTestFUN`. Must be one of "two.sided" (default), "greater" or "less".
#'
#' @details The element 'p.value' of the output is a `m x B` matrix whose entry i,j corresponds to `p_{(i)}(g_j.X)` with notation of the AoS 2020 paper cited below (section 4.5).
#'
#' @export
#'
#' @examples
#' m <- 50
#' n <- 45
#' X <- matrix(rnorm(m * n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#'
#' B <- 10
#' set.seed(123)
#' perm0 <- sanssouci:::get_randomized_p_values_two_sample(X, categ, B, rowWelchTests)
#'
#' # for this particular test 'get_randomized_p_values_two_sample' can be bypassed
#' set.seed(123)
#' null_groups <- replicate(B, sample(categ))
#' perm <- rowWelchTests(X, null_groups)
#' identical(perm0, perm$p.value)
#'
#' # Wilcoxon tests
#' set.seed(123)
#' perm0 <- sanssouci:::get_randomized_p_values_two_sample(X, categ, B, rowWilcoxonTests)
#' perm <- rowWilcoxonTests(X, null_groups)
#' identical(perm0, perm$p.value)
get_randomized_p_values_two_sample <- function(X, categ, B,
                                               rowTestFUN = rowWelchTests,
                                               alternative = c("two.sided", "less", "greater")) {
  m <- nrow(X)

  pval0 <- matrix(NA_real_, nrow = m, ncol = B)
  stat0 <- matrix(NA_real_, nrow = m, ncol = B)
  for (bb in 1:B) {
    categ_perm <- sample(categ)
    rwt <- rowTestFUN(X, categ = categ_perm, alternative = alternative)
    pval0[, bb] <- rwt$p.value
    stat0[, bb] <- rwt$statistic
  }

  return(pval0)
}

#' Get permutation  p-values
#'
#' Get a matrix of p-values under the null hypothesis obtained by sign-flipping (one-sample test).
#'
#' @param X A matrix of `m` variables (hypotheses) by `n` observations
#' @param B An integer value, the number of permutations to be performed
#' @param rowTestFUN a vectorized testing function (same I/O as [rowWelchTests])
#' @param seed a single value, interpreted as an integer or NULL. The seed used for the Random Number Generation. NULL by default.
#'
#' @details The element 'p.value' of the output is a `m x B` matrix whose entry i,j corresponds to `p_{(i)}(g_j.X)` with notation of the AoS 2020 paper cited below (section 4.5).
#'
#' @export
get_randomized_p_values_one_sample <- function(X, B = 100, rowTestFUN = rowTtestsOneSample, seed = NULL) {
  set.seed(seed)

  # Init
  n <- dim(X)[2]
  m <- dim(X)[1]

  # Initialize p-values
  pval0 <- matrix(0, nrow = m, ncol = B)
  # stat0 = matrix(0, nrow = m, ncol=B)

  for (bb in 1:B) {
    # X_flipped = t(t(X)*sample(x=c(-1,1), size=n, replace=TRUE))
    # rtonesample <- rowTestFUN(X_flipped)
    # pval0[,bb] <- rtonesample$p.value
    # stat0[,bb] <- rtonesample$statistic

    categ_permuted <- sample(x = c(-1, 1), size = n, replace = TRUE)
    ztest <- rowTestFUN(X, categ_permuted, alternative = "two.sided")
    pval0[, bb] <- ztest$p.value
    # stat0[,bb] <- ztest$statistic
  }

  pval0 <- apply(pval0, 2, sort)

  return(pval0)
}

#' Get a vector of pivotal statistics associated
#'   to permutation p-values and to a reference family
#'
#' @param p0 A matrix with B columns. Each row is a vector of m null p-values
#' @param m The total number of tested hypotheses
#' @param t_inv  An inverse threshold function (same I/O as 't_inv_linear')
#' @param K An integer value in `[1,m]`, the number of elements in the reference family. Defaults to `m`
#'
#' @return A vector of length `B` pivotal statitics, whose j-th entry
#'   corresponds to `psi(g_j.X)` with notation of the AoS 2020 paper cited
#'   below (section 4.5)
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#' @export
#'
#' @importFrom matrixStats colMins
#' @examples
#'
#' m <- 50
#' n <- 45
#' X <- matrix(rnorm(m * n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' null_groups <- replicate(B, sample(categ))
#' p0 <- rowWelchTests(X, null_groups)$p.value
#' pivStat <- get_pivotal_stat(p0, m)
#' quantile(pivStat, 0.2)
#'
#' @export
get_pivotal_stat <- function(p0, m,
                             t_inv = t_inv_linear,
                             K = nrow(p0)) {
  K <- force(K)
  stopifnot(m >= nrow(p0))
  stopifnot(K <= m)

  ## Step 2: order each column
  # p0 <- colSort(p0)
  p0s <- partialColSort(p0, K)
  B <- ncol(p0s)

  ## Step 3: apply template function
  tkInv <- matrix(nrow = K, ncol = B)
  for (kk in 1:K) {
    tkInv[kk, ] <- t_inv(p0s[kk, ], kk, m)
  }

  ## Step 4: report min for each column
  matrixStats::colMins(tkInv)
}


t_inv_linear <- function(y, k, m) {
  y * m / k
}

t_inv_beta <- function(y, k, m) {
  pbeta(y, k, m + 1 - k)
}

#' Get a vector of pivotal statistics associated
#'   to permutation p-values and to a reference family
#'
#' @param template matrix (m,B), template used for the calibration
#' @param pval0 matrix (m,B) of permuted pvalues
#' @param k_max An integer value in `[1,m]`, the number of elements in the reference family.
#'
#' @return A vector of length `B` pivotal statitics, whose j-th entry
#'   corresponds to `psi(g_j.X)` with notation of the AoS 2020 paper cited
#'   below (section 4.5)
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#' @export
estimate_jer <- function(template, pval0, k_max) {
  B <- dim(pval0)[2]
  m <- dim(pval0)[1]
  id_ranks <- matrix(rep(0:(m - 1), B), ncol = B, byrow = FALSE)

  cutoffs <- matrix(findInterval(pval0, template, left.open = TRUE), nrow = m)
  # print(cutoffs[k_max,])
  signs <- sign(id_ranks - cutoffs)
  sgn_trunc <- signs[1:(k_max),]

  JER <- mean(apply(sgn_trunc >= 0, 2, any))
  return(JER)
}

#' For a given risk level, calibrate the method on learned templates by
#'   dichotomy. This is equivalent to calibrating using pivotal stats but does
#'   not require the availability of a closed form inverse template.
#'
#' @param alpha A float, confidence level in \[0, 1\]
#' @param learned_templates A matrix with m rows and B columns. Learned templates for B permutations and m features
#' @param pval0 A matrix with m rows and B columns of permuted p-values.
#' @param k_max An integer value, the template size.
#' @param min_dist An integer, minimum distance to stop iterating dichotomy. Defaults to `m`
#'
#' @return A vector of length k_max. Threshold family chosen by calibration
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#' @export
#'
#' @importFrom matrixStats colMins
dichotomy <- function(alpha, learned_templates, pval0, k_max, min_dist = 3) {
  B <- dim(learned_templates)[2]
  m <- dim(learned_templates)[1]
  low <- 1
  high <- B

  if (estimate_jer(learned_templates[,high], pval0, k_max) <= alpha) {
    warning("All templates control the JER: choice may be conservative")
    return(learned_templates[1:k_max,high])
  }

  if (estimate_jer(learned_templates[, low], pval0, k_max) >= alpha) {
    warning("No suitable template found; Simes is used instead")
    piv_stat <- get_pivotal_stat(pval0, m = m, K = k_max)
    lambda <- stats::quantile(piv_stat, alpha, type = 1)
    thr <- t_linear(lambda, 1:k_max, m)
    return(thr)
  }


  while (high - low > min_dist) {
    mid <- as.integer((high + low) / 2)
    lw <- estimate_jer(learned_templates[, low], pval0, k_max) - alpha
    md <- estimate_jer(learned_templates[, mid], pval0, k_max) - alpha
    if (md == 0) {
      return(mid)
    } else if (lw * md < 0) {
      high <- mid
    } else {
      low <- mid
    }
  }
  return(learned_templates[1:k_max,low])
}

#' For a given risk level, calibrate the method on learned templates by
#'   dichotomy. This is equivalent to calibrating using pivotal stats but does
#'   not require the availability of a closed form inverse template.
#'
#' @param X A matrix (m, n)
#' @param label A numeric vector of \eqn{n} values in \eqn{0, 1}, the groups of
#'   observations on which to perform two-sample tests
#' @param B integer, the number of permutation.
#' @param row_test_fun A (vectorized) test function. Defaults to [rowWelchTests].
#' @param alternative A character string specifying the alternative hypothesis, to be passed to `rowTestFUN`. Must be one of "two.sided" (default), "greater" or "less".
#'
#' @return A matrix (m, B) of a learned template
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#'
#' @importFrom matrixStats colMins
#' @examples
#'
#' m <- 50
#' n <- 45
#' X <- matrix(rnorm(m * n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' null_groups <- replicate(B, sample(categ))
#' p0 <- rowWelchTests(X, null_groups)$p.value
#' pivStat <- get_pivotal_stat(p0, m)
#' quantile(pivStat, 0.2)
#'
get_data_driven_template <- function(X, label, B,
                                     row_test_fun = rowWelchTests,
                                     alternative = c("two.sided", "less", "greater")) {
  learned_template_ <- get_randomized_p_values_two_sample(
    X = X, categ = label, B = B,
    rowTestFUN = row_test_fun,
    alternative = alternative
  )
  learned_template_ <- t(apply(learned_template_, 1, sort))
  learned_template <- apply(learned_template_, 2, sort)
  return(learned_template)
}



#' Perform JER calibration
#'
#' @param X A matrix of `m` variables (hypotheses) by `n` observations
#' @param label A numeric vector of \eqn{n} values in \eqn{0, 1}, the groups of
#'   observations on which to perform two-sample tests
#' @param B An integer value, the number of permutations to be performed
#' @param row_test_fun A (vectorized) test function. Defaults to [rowWelchTests].
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#' @param template A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes", or "Oracle". "Linear" and "Simes" families are
#'   identical.
#' @param k_max A numeric value in `[0,m]`, the length of the template
#'
#' @return A list with elements
#' - thr: A numeric vector of length K, such that the estimated probability that
#' there exists an index k between 1 and K such that the k-th maximum of the
#' test statistics of is greater than `thr[k]`, is less than alpha
#' - piv_stat: A vector of `B` pivotal statitics
#' - lambda: A numeric value, the result of the calibration
#'
#' @details 'calibrate0' performs single step calibration,
#'   whereas 'calibrate' performs step-down calibration. Hence
#'   the output of 'calibrate(..., max_steps_down = 0)' and
#'   'calibrate0(...)' should be identical
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. *Annals of
#'   Statistics*, 48(3), 1281-1303.
#'
#' @export
calibrate_jer <- function(X, label, alpha, B = 100, row_test_fun = rowWelchTests, template = "Simes", k_max = nrow(X)) {
  if (is.null(label)) {
    perm_p_values <- get_randomized_p_values_one_sample(X, B = B)
  } else {
    categCheck(label, ncol(X))
    perm <- get_randomized_p_values_two_sample(X, label, B, rowTestFUN = row_test_fun, alternative = "two.sided")
    perm_p_value <- perm$p.value
  }

  if (any(template == c("Simes", "Linear", "Beta"))) {
    calibration <- calibrate0(p0 = perm_p_value, m = nrow(X), alpha = alpha, family = template, K = k_max)
    return(calibration$thr[calibration$lambda]) # ?
  } else if (is.numeric(template)) {
    lambda <- dichotomy(alpha = alpha, learned_templates = template, pval0 = perm_p_value, k_max = k_max)
    return(template[lambda])
  }
}


#' Perform JER calibration with the Simes template for one sample.
#'
#' @param X A matrix of `m` variables (hypotheses) by `n` observations
#' @param B An integer value, the number of permutations to be performed
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#'
#' @return A list with elements
#' - thr: A numeric vector of length K, such that the estimated probability that
#' there exists an index k between 1 and K such that the k-th maximum of the
#' test statistics of is greater than `thr[k]`, is less than alpha
#' - piv_stat: A vector of `B` pivotal statitics
#' - lambda: A numeric value, the result of the calibration
#'
#' @details calibrate_simes_one_sample run calibrate_jer for one sample.
#' Label is NULL, the template is fixed to "Simes" and k_max = nrow(X)
calibrate_simes_one_sample <- function(X, alpha, B = 100) {
  # ...
  return(calibrate_jer(X = X, label = NULL, alpha = alpha, B = B, template = "Simes", k_max = nrow(X)))
}

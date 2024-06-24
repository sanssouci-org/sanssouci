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
#' sim <- gaussianSamples(m = m, rho = 0.3, n = 45, 
#'                        pi0 = 0.8, SNR = 3, prob = 0.5)
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
#' X <- expr_ALL; rm(expr_ALL)
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
#' minTP(p, thr)  ## lower bound on true positives
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
#' FDP_bound <- sanssouci:::curveMaxFP(sort(p), thr)/seq(along = p)
#' plot(head(FDP_bound, 300), t = 's', 
#'   xlab = "Number of features",
#'   ylab = "Upper bound on False Discovery Proportion")
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
                      piv_stat0 = piv_stat0)
    thr <- cal$thr[1]   ## (1-)FWER threshold
    R1 <- integer(0)
    if (!is.null(p)) {
        R1 <- which(p <= thr)
    }
    
    ## force 'convergence' if nothing to gain) or nothing left to be rejected
    converged <- (length(R1) == 0L)  | (length(R1) == m)
    
    while (!converged && step < max_steps_down) {
        step <- step + 1
        
        p1 <- p0[-R1, ]
        thr <- cal$thr
        lambda <- cal$lambda
        cal <- calibrate(p1, m, alpha, 
                         family = family, 
                         K = K,
                         piv_stat0 = NULL) ## force piv stat calc
        R1_new <- which(p < cal$thr[1])
        
        noNewRejection <- all(R1_new %in% R1)          ## convergence reached?
        if (noNewRejection) {
            if (!identical(R1_new, R1)) {              ## can this happen?? covr check
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
    lambda <- stats::quantile(pivStat, alpha, type = 1) # 'type=1' corresponds to inverse empirical cdf as per (20) in BNR 2020
    thr <- t_(lambda, 1:K, m)
    res <- list(thr = thr,
                piv_stat = pivStat,
                lambda = lambda)
    return(res)
}

#' Get permutation statistics and p-values
#' 
#' Get a matrix of statistics and p-values under the null hypothesis obtained by repeated permutation of class labels (two-sample tests)
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
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' 
#' B <- 10
#' set.seed(123)
#' perm0 <- sanssouci:::get_perm(X, categ, B, rowWelchTests)
#' 
#' # for this particular test 'get_perm' can be bypassed
#' set.seed(123)
#' null_groups <- replicate(B, sample(categ))
#' perm <- rowWelchTests(X, null_groups)
#' identical(perm0$p.value, perm$p.value)
#' identical(perm0$statistic, perm$statistic)

#' # Wilcoxon tests
#' set.seed(123)
#' perm0 <- sanssouci:::get_perm(X, categ, B, rowWilcoxonTests)
#' perm <- rowWilcoxonTests(X, null_groups)
#' identical(perm0$p.value, perm$p.value)
#' identical(perm0$statistic, perm$statistic)
get_perm <- function(X, categ, B, 
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
    
    list(p.value = pval0,
         statistic = stat0)
}

#' Get a vector of pivotal statistics associated
#'   to permutation p-values and to a reference family
#'
#' @param p0 A matrix with B rows. Each row is a vector of null p-values
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
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' null_groups <- replicate(B, sample(categ))
#' p0 <- rowWelchTests(X, null_groups)$p.value
#' pivStat <- get_pivotal_stat(p0, m)
#' quantile(pivStat, 0.2)
#' 
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
    y * m / k;
}

t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}



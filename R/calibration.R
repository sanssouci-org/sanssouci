#' Get permutation statistics and p-values
#' 
#' Get a matrix of statistics and p-values under the null hypothesis obtained by repeated permutation of class labels (two-sample tests)
#' 
#' @param X A matrix of `m` variables (hypotheses) by `n` observations
#' @param categ An numeric vector of `n` values in `0,1`
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#' @param B An integer value, the number of permutations to be performed
#' @param alpha A numeric value in `[0,1]`, the target risk
#' @param rowTestFUN a vectorized testing function (same I/O as [rowWelchTests])
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
#' perm0 <- sansSouci:::get_perm(X, categ, B, rowWelchTests)
#' 
#' # for this particular test 'get_perm' can be bypassed
#' set.seed(123)
#' null_groups <- replicate(B, sample(categ))
#' perm <- rowWelchTests(X, null_groups)
#' identical(perm0$p.value, perm$p.value)
#' identical(perm0$statistic, perm$statistic)

#' # Wilcoxon tests
#' set.seed(123)
#' perm0 <- sansSouci:::get_perm(X, categ, B, rowWilcoxonTests)
#' perm <- rowWilcoxonTests(X, null_groups)
#' identical(perm0$p.value, perm$p.value)
#' identical(perm0$statistic, perm$statistic)
get_perm <- function(X, categ, B, 
                       rowTestFUN = rowWelchTests) {
    m <- nrow(X)
    
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) 
    stat0 <- matrix(NA_real_, nrow = m, ncol = B) 
    for (bb in 1:B) {
        categ_perm <- sample(categ)
        rwt <- rowTestFUN(X, categ = categ_perm)
        pval0[, bb] <- rwt$p.value
        stat0[, bb] <- rwt$statistic
    }
    
    list(p.value = pval0,
         statistic = stat0)
}

#' Get a vector of pivotal statistics associated
#'   to permutation p-values and to a reference family
#'
#' @param alpha A numeric value in `[0,1]`, the target risk
#' @param family A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes". "Linear" and "Simes" families are
#'   identical.
#' @param p0 A `m x B` matrix. The j-th ow corresponds to a permutation
#'   of the input `categ`: for each hypothesis i, `p0[i,j]` is the
#'   p-value of the test of the i-th null hypothesis on the permuted categories
#' @param t_inv  An inverse threshold function (same I/O as [t_inv_linear])
#' @param m The total numer of tested hypotheses. Defaults to `p0`
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
#' pivStat <- get_pivotal_stat(p0)
#' quantile(pivStat, 0.2)
#' 
#' 
#' @export
get_pivotal_stat <- function(p0,
                             t_inv = t_inv_linear,
                             m = nrow(p0),
                             K = m) {
    m <- force(m)
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


#' Perform JER calibration from randomization p-values
#' 
#' @inheritParams get_pivotal_stat
#' @param family A character vector, the name of the reference family. Should be
#'   either "Simes" (aka "Linear") or "Beta"
#' @param piv_stat0 Don't use! Should be removed soon...
#' @param max_steps_down A numeric value, the maximum number of steps down to
#'   perform. Defaults to 10 (but the algorithm generally converges in 1 or 2
#'   steps).
#' 
#' @return A list with elements
#' - thr: A numeric vector of length K, such that the estimated probability that there exists an index k between 1 and K such that the k-th maximum of the test statistics of is greater than `thr[k]`, is less than alpha
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
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 100
#' p <- rowWelchTests(X, categ)$p.value
#' 
#' null_groups <- replicate(B, sample(categ))
#' p0 <- rowWelchTests(X, null_groups)$p.value
#' calib0 <- calibrate0(p0, alpha = 0.1) # single step
#' calib <- calibrate(p0, alpha = 0.1)
#' calib$lambda >= calib0$lambda # probably very close here (null data)
#' 
#' maxFP(p, calib$thr)
#' 
#' \dontrun{
#' # Gene expression
#' data(expr_ALL, package = "sansSouci.data")
#' X <- expr_ALL; rm(expr_ALL)
#' categ <- ifelse(colnames(X) == "BCR/ABL", 1, 0) # map to 0/1
#' 
#' null_groups <- replicate(500, sample(categ))
#' perm <- rowWelchTests(X, null_groups)
#' p0 <- perm$p.value
#' 
#' alpha <- 0.1
#' calib_L <- calibrate(p0, alpha, family = "Linear")
#' calib_B <- calibrate(p0, alpha, family = "Beta", K = 100)
#' p <- rowWelchTests(X, categ)$p.value
#' 
#' ## post hoc bounds (these are functions!)
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
#' FDP_bound <- sansSouci:::curveMaxFP(sort(p), thr)/seq(along = p)
#' plot(head(FDP_bound, 300), t = 's', 
#'   xlab = "Number of features",
#'   ylab = "Upper bound on False Discovery Proportion")
#' }
#' 
#' @export
calibrate <- function(p0, alpha, 
                      family = c("Linear", "Beta", "Simes"), 
                      m = nrow(p0),
                      K = m,
                      p = NULL, 
                      max_steps_down = 10L,
                      piv_stat0 = NULL) {
    step <- 0
    cal <- calibrate0(p0, alpha, 
                      family = family, 
                      m = m,
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
        cal <- calibrate(p1, alpha, 
                         family = family, 
                         m = m,
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
calibrate0 <- function(p0, alpha, 
                       family = c("Linear", "Beta", "Simes"), 
                       m = nrow(p0),
                       K = m,
                       piv_stat0 = NULL) {
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
        pivStat <- get_pivotal_stat(p0, t_inv, m, min(nrow(p0), K))
    }
    lambda <- stats::quantile(pivStat, alpha, type = 1)
    thr <- t_(lambda, 1:K, m)
    res <- list(thr = thr,
                piv_stat = pivStat,
                lambda = lambda)
    return(res)
}

t_inv_linear <- function(y, k, m) {
    y * m / k;
}

t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}


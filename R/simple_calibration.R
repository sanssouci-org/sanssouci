t_inv_linear <- function(y, k, m) {
    y * m / k;
}

t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}


#' @title Low-level calibration functions
#' 
#' @description Get a matrix of p-values under the null hypothesis obtained 
#' by repeated permutation of class labels
#' 
#' @rdname calibrate_low-level
#' @inheritParams calibrate
#' 
#' @return A \eqn{m \times B} matrix whose entry i,j corresponds to \eqn{p_{(i)}(g_j.X)} with notation of the AoS 2020 paper cited below (section 4.5) 
#' @export
#' 
#' @examples
#' 
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' p0 <- sansSouci:::get_perm_p(X, categ, B)
get_perm_p <- function(X, categ, B, 
                       rowTestFUN = rowWelchTests) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values for B permutations of the class assignments
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) ## parametric p-values
    for (bb in 1:B) {
#        categ_perm <- sample(categ, length(categ))
        categ_perm <- sample(categ)
        rwt <- rowTestFUN(X, categ = categ_perm)
        pval0[, bb] <- rwt$p.value
    }
    
    ## Step 2: sort each column
    return(colSort(pval0))
}

#' @description get_pivotal_stat: Get a vector of pivotal statistics associated to permutation p-values and to a reference family
#' 
#' @rdname calibrate_low-level
#' @inheritParams calibrate
#' @param p0 A \eqn{m \times B} matrix. The j-th ow corresponds to a permutation
#'   of the input \code{categ}: for each hypothesis i, \code{p0[i,j]} is the
#'   p-value of the test of the i-th null hypothesis on the permuted categories
#' @param t_inv  An inverse threshold function (same I/O as \code{t_inv_linear})
#'   
#' @return A vector of length \eqn{B} pivotal statitics, whose j-th entry
#'   corresponds to \eqn{\psi(g_j.X)} with notation of the AoS 2020 paper cited
#'   below (section 4.5)
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families. Annals of
#'   Statistics, 48(3), 1281-1303.
#'   
#' @export
#' 
#' @importFrom matrixStats colMins
#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' p0 <- sansSouci:::get_perm_p(X, categ, B)
#' pivStat <- sansSouci:::get_pivotal_stat(p0)
#' 
get_pivotal_stat <- function(p0,
                             t_inv = t_inv_linear,
                             K = nrow(p0)) {
    m <- nrow(p0)
    B <- ncol(p0)

    ## Step 3: apply template function
    tkInv <- matrix(nrow = K, ncol = B)
    for (kk in 1:K) {
        tkInv[kk, ] <- t_inv(p0[kk, ], kk, m)
    }
    
    ## Step 4: report min for each column
    matrixStats::colMins(tkInv)
}

#' @description Get a vector of calibrated thresholds associated to a pivotal statistics and a reference family
#' 
#' @rdname calibrate_low-level
#' @inheritParams calibrate
#' 
#' @return A JER controlling family, that is, a numeric vector of length
#'   \code{K}, such that the estimated probability that there exists an index
#'   \eqn{k} between 1 and \eqn{K} such that the \eqn{k}-th maximum of the test
#'   statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}
#' 
#' @export
#' 
#' @examples
#' 
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 100
#' p <- rowWelchTests(X, categ)$p.value
#' p0 <- sansSouci:::get_perm_p(X, categ, B)
#' calib <- get_calibrated_thresholds(p0, alpha = 0.1, family = "Linear")
#' #Vbar <- get_max_FP(p, calib$thr)
#' #Vbar(1:m)
#' 
#' # Gene expression
#' data(expr_ALL, package = "sansSouci.data")
#' categ <- ifelse(colnames(expr_ALL) == "BCR/ABL", 1, 0) # map to 0/1
#' 
#' p0 <- sansSouci:::get_perm_p(expr_ALL, categ, B = 100)
#' alpha <- 0.1
#' calib_L <- get_calibrated_thresholds(p0, alpha, family = "Linear")
#' calib_B <- get_calibrated_thresholds(p0, alpha, family = "Beta")
#' p <- rowWelchTests(expr_ALL, categ)$p.value
#' 
#' ## post hoc bounds (these are functions!)
#' #FP <- get_max_FP(p, calib_L$thr)  ## upper bound on false positives
#' #FDP <- get_max_FDP(p, calib_L$thr) ## upper bound on FDP
#' 
#' ## example of user selection: rejections of BH(0.05) procedure

#' adjp <- p.adjust(p, method = "BH") 
#' sel <- which(adjp < 0.05)
#' length(sel)
#' 
#' #FP(sel)
#' #FDP(sel)
#' 
#' # confidence envelope on the FDP
#' #FDP_env <- FDP(seq_along(p), envelope = TRUE)
#' #plot(head(FDP_env, 300), t = 's', 
#' #  xlab = "Number of features",
#' #  ylab = "Upper bound on FDP")
#' 
#' 
get_calibrated_thresholds <- function(p0, alpha, 
                                      family = c("Linear", "Beta", "Simes"), 
                                      K = nrow(p0)) {
    family <- match.arg(family)
    if (family %in% c("Linear", "Simes")) {
        t_inv <- t_inv_linear
        t_ <- t_linear
    } else if (family == "Beta") {
        t_inv <- t_inv_beta
        t_ <- t_beta
    }
    pivStat <- get_pivotal_stat(p0, t_inv, K)
    lambda <- stats::quantile(pivStat, alpha, type = 1)
    thr <- t_(lambda, 1:K, nrow(p0))
    res <- list(thr = thr,
                pivStat = pivStat,
                lambda = lambda)
    return(res)
}


#' Calibrate 
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param categ An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#' @param B An integer value, the number of permutations to be performed
#' @param alpha A numeric value in [0,1], the target risk
#' @param rowTestFUN a vectorized testing function (same I/O as \code{rowWelchTests})
#' @param family A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes". "Linear" and "Simes" families are 
#'   identical.
#' @param K An integer value in \eqn{[1,m]}, the number of elements in the
#'   reference family. Defaults to \eqn{m}.
#' @return A list with elements \describe{
#'   \item{thr}{}
#'   \item{p0}{A \eqn{B \times m} matrix of p-values under the null hypothesis}
#'   \item{pivStat}{A vector of \eqn{B} pivotal statitics}
#'   \item{lambda}{A numeric value, the result of the calibration} 
#' }
#' @export
#' 
#' @examples
#'
#' # 1. Simulated data
#' set.seed(0xBEEF)
#' m <- 500
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pval <-  rowWelchTests(X, categ)$p.value
#' cal <- calibrate(X, categ, B = 1e2, alpha = 0.1, family = "Simes")
#' 
#' # Lower bound on number of true DEG in entire data set
#' #TP <- get_min_TP(pval, cal$thr)
#' #TP(1:m)
#' 
#' # 2. Gene expression
#' data(expr_ALL, package = "sansSouci.data")
#' X <- expr_ALL
#' m <- nrow(X)
#' categ <- ifelse(colnames(X) == "BCR/ABL", 1, 0) # map to 0/1
#' p <-  rowWelchTests(X, categ)$p.value
#' alpha <- 0.1
#' cal <- calibrate(X, categ, B = 1e2, alpha = 0.1, family = "Simes")
#' 
#' ## post hoc bounds (these are functions!)
#' #FP <- get_max_FP(p, cal$thr)  ## upper bound on false positives
#' #FDP <- get_max_FDP(p, cal$thr) ## upper bound on FDP
#' 
#' ## example of user selection: rejections of BH(0.05) procedure
#' adjp <- p.adjust(p, method = "BH") 
#' sel <- which(adjp < 0.05)
#' length(sel)
#' 
#' #FP(sel)
#' #FDP(sel)
#' 
#' # confidence envelope on the FDP
#' #FDP_env <- FDP(seq_along(p), envelope = TRUE)
#' #plot(head(FDP_env, 300), t = 's', 
#' #  xlab = "Number of features",
#' #  ylab = "Upper bound on FDP")
#'   
#' # Lower bound on number of true DEG in entire data set
#' #TP <- get_min_TP(p, cal$thr)
#' #TP(1:m)
#' 
#' #FP <- get_max_FP(p, cal$thr)
#' #FP(1:m)
#' 
#' # Compare to Simes (without calibration)
#' #FDP_Simes <- get_max_FDP(p, alpha*1:m/m)
#' #bound <- FDP_Simes(1:m, envelope = TRUE)
#' #lines(head(bound, 300), t = "s", lty = 2)
#' 
calibrate <- function(X, categ, B, alpha,
                      rowTestFUN = rowWelchTests, 
                      family = c("Linear", "Beta", "Simes"), 
                      K = nrow(X)) {
    p0 <- get_perm_p(X, categ, B = B, rowTestFUN = rowTestFUN)
    calib <- get_calibrated_thresholds(p0, alpha, 
                                     family = family, K = K)
    res <- list(thr = calib$thr,
                p0  = p0,
                pivStat = calib$pivStat,
                lambda = calib$lambda)
    return(res)
}

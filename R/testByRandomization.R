#' Randomization-based testing
#'
#' Randomization-based testing using permutation or sign-flipping
#'
#' @param X a matrix of \code{m} variables by \code{n} observations
#'
#' @param B A numeric value, the number of permutations to be performed
#'
#' @param cls A vector of length \code{n} class labels in \code{0,1} for flavor
#'   "perm". Defaults to colnames(X).
#'
#' @param p.value A boolean value: should randomization \eqn{p}-values be
#'   calculated and returned? Defaults to @TRUE
#'
#' @param seed An integer (or NULL) value used as a seed for random number
#'   generation. If \code{NULL}, no seed is specified
#'
#' @details The type of randomization is determined by the parameter \code{cls}.
#'   If \code{cls} does not containt two distinct values (or is \code{NULL}), a
#'   one-sample test is performed using randomization (flavor "flip"). If it
#'   contains two distinct values, a two-sample test is perfomed using
#'   permutations (flavor "perm").
#'
#'   For permutation, we test the null hypothesis: "both groups have the same
#'   mean" against the one-sided alternative that the mean is larger in the
#'   second group. The test is Welch's two-sample test for unequal variances.
#'   Permuted test statistics are calculated by B permutations of the group
#'   labels. Corresponding observed and permuted p-values are calculated as the
#'   proportion of permutations (including the identity) for which the permuted
#'   test statistic is larger than the observed test statistic. The
#'   corresponding parametric p-values are returned as elements \code{param.p}
#'   (for the original data) and \code{param.p0} (for the permutations).
#'
#'   For sign-flipping, we test the null hypothesis: "the mean is 0" against the
#'   two-sided alternative that the mean is larger than 0. We use the (rescaled)
#'   empirical mean of the observations as a test statistic. Sign-flipped test
#'   statistics are calculated by flipping the sign of each observation with
#'   probability 1/2.
#'
#' @references Ge, Y., Dudoit, S. and Speed, T.P., 2003. Resampling-based
#'   multiple testing for microarray data analysis. _Test_, 12(1), pp.1-77.
#'
#' @return a list with elements: \describe{
#'
#'   \item{T}{A vector of \eqn{m} test statistics}
#'
#'   \item{T0}{A \eqn{m \times B} matrix of randomized test statistics}
#'
#'   \item{p}{A vector of \eqn{m} \eqn{p}-values (only if \code{p.value} is
#'   \code{TRUE} )}
#'
#'   \item{p0}{A \eqn{m \times B} matrix of randomization \eqn{p}-values (only
#'   if \code{p.value} is \code{TRUE} )}
#'
#'   \item{flavor}{A character value, the type of randomization performed:
#'   "perm" for permutation-based randomization in two-sample tests, and "flip"
#'   for sign-flipping-based randomization in one sample tests. See Details.}
#'
#'   \item{param.p}{A vector of \eqn{m} parametric \eqn{p}-values (only for
#'   flavor "perm")}
#'
#'   \item{param.p0}{A \eqn{m \times B} matrix of parametric \eqn{p}-values on
#'   permuted data (only for flavor "perm" )}
#'
#'   \item{df}{A vector of \eqn{m} degrees of freedom for the observed
#'   statistics (only for flavor "perm")}
#'
#'   \item{df0}{A \eqn{m \times B} matrix of degrees of freedom on permuted data
#'   (only for flavor "perm" )}}
#'
#'
#' @examples
#'
#' p <- 510
#' n <- 380
#' B <- 1e3
#' mat <- matrix(rnorm(p*n), ncol=n)
#' cls <- rep(c(0, 1), times=c(n/2, n-n/2))
#' resPerm <- testByRandomization(X=mat, cls=cls, B=B) ## permutation
#' resFlip <- testByRandomization(X=mat, B=B)          ## sign-flipping
#'
#' # empirical coverage of Simes thresholds
#' alpha <- 0.2
#' thr <- SimesThresholdFamily(p)(alpha)
#' sansSouci:::empiricalCoverage(thr, resPerm$T0) ## Welch, not Gaussian
#' sansSouci:::empiricalCoverage(thr, qnorm(1-resPerm$p0))
#' sansSouci:::empiricalCoverage(thr, qnorm(1-resPerm$param.p0))
#'
#' sansSouci:::empiricalCoverage(thr, resFlip$T0)
#' sansSouci:::empiricalCoverage(thr, qnorm(1-resFlip$p0))
#'
#' # test statistics null distribution
#' hist(resPerm$p)
#' hist(resFlip$p)
#'
#' alpha <- 0.05
#' resP <- jointFWERControl(resPerm$T0, refFamily="Simes", alpha=alpha, stat=resPerm$T)
#' resF <- jointFWERControl(resFlip$T0, refFamily="Simes", alpha=alpha, stat=resFlip$T)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
testByRandomization <- function(X, B, cls = colnames(X), 
                                p.value = TRUE, seed = NULL){
    ## sanity checks
    n <- ncol(X)
    luc <- length(unique(cls))
    if (luc <= 1) {  
        # no classes or a single class given: assuming sign flipping 
        flavor <- "flip"
    } else {
        if (length(cls) != n) {
            stop("The number of columns of argument 'X' should match the length of argument 'cls'")
        }
        if (luc == 2) {
            flavor <- "perm"
            tbl <- table(cls)
            if ( !all(names(tbl) == c("0", "1"))) {  # note that numeric values are allowed in cls as they are converted into character by 'table'...
                stop("Argument 'cls' should contain (only) 0:s and 1:s")
            }
            if (min(tbl) < 3) {
                stop("Argument 'cls' should contain at least 3 elements of each sample")
            }
        } else if (luc > 2) {
            stop("Tests for more than 2 classes not implemented yet")
        }
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    m <- nrow(X)
    
    if (flavor == "perm") {
        ## TODO: (cf issue #3)
        ## * other statistics ? (difference in empirical means, Mann-Whitney)
        ## * one-sided tests ?
        
        ## observed test statistics
        rwt <- rowWelchTests(X, categ = cls)
        T_obs <- rwt$statistic
        p_obs <- rwt$p.value  ## parametric p-value
        df_obs <- rwt$parameter  ## degrees of freedom of the T statistics
        ## test statistics under H0
        T <- matrix(nrow = m, ncol = B)
        pp <- matrix(nrow = m, ncol = B) ## parametric p-value
        df <- matrix(nrow = m, ncol = B) 
        for (bb in 1:B) {
            cls_perm <- sample(cls, length(cls))
            rwt <- rowWelchTests(X, categ = cls_perm)
            T[, bb] <- rwt$statistic
            pp[, bb] <- rwt$p.value
            df[, bb] <- rwt$parameter
        }
        res <- list(T = T_obs, T0 = T, 
                    flavor = flavor,
                    param.p = p_obs, param.p0 = pp,
                    df = df_obs, df0 = df)
    } else if (flavor == "flip") {
        ## observed test statistics
        T_obs <- rowSums(X)/sqrt(n)
        ## test statistics under H0
        T <- testBySignFlipping(X, B)
        res <- list(T = T_obs, T0 = T, flavor = flavor)
    }

    if (p.value) {
        ## get m x (B+1) matrix of pvalues under the null (+ original)
        ## by sorting null test statistics as proposed by Ge et al (2003)
        TT <- cbind(T, T_obs)
        pB <- rowRanks(-abs(TT)) / (B+1)
        
        res$p <- pB[, B+1]
        res$p0 <- pB[, -(B+1), drop = FALSE]
    }
    return(res)
}


testBySignFlippingR <- function(X, B) {
    m <- nrow(X)
    n <- ncol(X)
    
    T <- matrix(nrow = m, ncol = B)
    for (bb in 1:B) {
        eps <- rbinom(n, 1, prob = 0.5)*2 - 1  ## signs
        eX <- sweep(X, MARGIN = 2, STATS = eps, FUN = `*`)
        Tb <- rowSums(eX)/sqrt(n)
        T[, bb] <- Tb
    }
    T
}

testByPermutationR <- function(X, cls, B) {
    m <- nrow(X)
    n <- ncol(X)
    stopifnot(n == length(cls))
    
    T <- matrix(nrow = m, ncol = B)
    for (bb in 1:B) {
        cls_perm <- sample(cls, length(cls))
        Tb <- rowWelchTests(X, categ = cls_perm)$statistic
        T[, bb] <- Tb
    }
    T
}
#' Randomization-based testing for equal means in two groups
#' 
#' @param X a matrix of \code{m} variables by \code{n} observations
#'   
#' @param B A numeric value, the number of permutations to be performed
#'   
#' @param flavor A character value, the type of randomization to be 
#'   performed. Should either be "perm" for two-sample permutation or "flip" for
#'   sign flipping
#'   
#' @param cls A vector of length \code{n} class labels in \code{0,1} for flavor
#'   "perm". Defaults to NULL
#'   
#' @param p.value A boolean value: should \eqn{p}-values be calculated and 
#'   returned? Defaults to @TRUE
#'   
#' @param seed An integer (or NULL) value used as a seed for random number 
#'   generation. If \code{NULL}, no seed is specified
#'   
#' @details The test statistics is Welch's two-sample test for unequal variances
#'   
#' @references Ge, Y., Dudoit, S. and Speed, T.P., 2003. Resampling-based 
#'   multiple testing for microarray data analysis. _Test_, 12(1), pp.1-77.
#'   
#' @return a list with two or four elements: \describe{ \item{T}{A vector of
#'   \eqn{m} test statistics} \item{T0}{A \eqn{m \times B} matrix of permutation
#'   test statistics} \item{p}{A vector of \eqn{m} \eqn{p}-values (only if
#'   \code{p.value} is @TRUE)} \item{p0}{A \eqn{m \times B} matrix of
#'   permutation \eqn{p}-values (ony if \code{p.value} is @TRUE)} }
#'   
#' @examples
#' 
#' p <- 513
#' n <- 38
#' mat <- matrix(rnorm(p*n), ncol=n)
#' cls <- rep(c(0, 1), times=c(27, n-27))
#' resPerm <- testByRandomization(X=mat, flavor="perm", cls=cls, B=1000)
#' resFlip <- testByRandomization(X=mat, flavor="flip", B=1000)
#' 
#' # show test statistics null distribution
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
testByRandomization <- function(X, B, flavor=c("perm", "flip"), cls=NULL, p.value=TRUE, seed=NULL){
    ## sanity checks
    n <- ncol(X)
    flavor <- match.arg(flavor)
    if (flavor == "perm") {
        if (length(cls) != n) {
            stop("The number of columns of argument 'X' should match the length of argument 'cls'")
        }
        if (!all(sort(unique(cls)) == c(0, 1))) {
            stop("Argument 'cls' should contain (only) 0:s and 1:s")
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
        T_obs <- rowWelchTests(X, categ = cls)$statistic
        ## test statistics under H0
        T <- matrix(nrow = m, ncol = B)
        for (bb in 1:B) {
            cls_perm <- sample(cls, length(cls))
            Tb <- rowWelchTests(X, categ = cls_perm)$statistic
            T[, bb] <- Tb
        }
    } else if (flavor == "flip") {
        ## observed test statistics
        T_obs <- rowSums(X)/sqrt(n)
        ## test statistics under H0
        T <- testBySignFlipping(X, B)
    }

    res <- list(T = T_obs, T0 = T)
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
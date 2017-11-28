#' Perform permutation testing for equal means in two groups
#' 
#' @param X a matrix of \code{m} variables by \code{n} observations
#' 
#' @param cls a vector of length \code{n} class labels in \code{0,1}
#'   
#' @param B a numeric value, the number of permutations to be performed
#' 
#' @param p.value a boolean value: should \eqn{p}-values be calculated and 
#'   returned? Defaults to @TRUE
#'   
#' @param seed an integer (or NULL) value used as a seed for random number 
#'   generation. If \code{NULL}, no seed is specified.
#'   
#' @details The test statistics is Welch's two-sample test for unequal variances
#'   
#' @references Ge, Y., Dudoit, S. and Speed, T.P., 2003. Resampling-based
#'   multiple testing for microarray data analysis. _Test_, 12(1), pp.1-77.
#'   
#' @return a list with two or four elements: \describe{
#' \item{T}{A vector of \eqn{m} test statistics}
#' \item{T0}{A \eqn{m \times B} matrix of permutation test statistics}
#' \item{p}{A vector of \eqn{m} \eqn{p}-values (only if \code{p.value} is @TRUE)}
#' \item{p0}{A \eqn{m \times B} matrix of permutation \eqn{p}-values (ony if \code{p.value} is @TRUE)}
#' }
#' 
#' @examples
#' 
#' p <- 513
#' n <- 38
#' mat <- matrix(rnorm(p*n), ncol=n)
#' cls <- rep(c(0, 1), times=c(27, n-27))
#' resPerm <- testByTwoSamplePermutation(X=mat, cls=cls, B=1000)
#' 
#' alpha <- 0.05
#' res <- jointFWERControl(resPerm$T0, refFamily="Simes", alpha=alpha, stat=-resPerm$T)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
testByTwoSamplePermutation <- function(X, cls, B, p.value=TRUE, seed=NULL){
    ## sanity checks
    n <- ncol(X)
    if (length(cls) != n) {
        stop("The number of columns of argument 'X' should match the length of argument 'cls'")
    }
    if (!all(sort(unique(cls)) == c(0, 1))) {
        stop("Argument 'cls' should contain (only) 0:s and 1:s")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    m <- nrow(X)
    T <- matrix(nrow=m, ncol=B+1)
    
    ## TODO: (cf issue #3)
    ## * other statistics ? (difference in empirical means, Mann-Whitney)
    ## * one-sided tests ?
    
    ## get original test statistics
    T_obs <- rowWelchTests(X, categ = cls)$statistic
    T[, B+1] <- T_obs
    ## get permutation test statistic
    for (bb in 1:B) {
        cls_perm <- sample(cls, length(cls))
        Tb <- rowWelchTests(X, categ = cls_perm)$statistic
        T[, bb] <- Tb
    }
    ## Compute the observed test stat
    
    res <- list(T = T_obs, T0 = T[, -(B+1), drop = FALSE])
    if (p.value) {
        ## get m x (B+1) matrix of pvalues under the null (+ original)
        ## by sorting null test statistics as proposed by Ge et al (2003)
        pB <- rowRanks(-abs(T)) / (B+1)

        res$p <- pB[, B+1]
        res$p0 <- pB[, -(B+1), drop = FALSE]
    }
    return(res)
}


#' Test by sign flipping in a location model
#' 
#' @param X a matrix of \code{m} variables by \code{n} observations
#' 
#' @param B a numeric value, the number of flips to be performed
#' 
#' @param p.value a boolean value: should \eqn{p}-values be calculated and 
#'   returned? Defaults to @TRUE
#'   
#' @param seed an integer (or NULL) value used as a seed for random number 
#'   generation. If \code{NULL}, no seed is specified.
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
#' resFlip <- testBySignFlipping(X=mat, B=1000, seed=123)
#' 
#' alpha <- 0.05
#' res <- jointFWERControl(resFlip$T0, refFamily="Simes", alpha=alpha, stat=resFlip$T)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
testBySignFlipping <- function(X, B, p.value=TRUE, seed=NULL){
    ## sanity checks
    n <- ncol(X)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    m <- nrow(X)
    T <- matrix(nrow = m, ncol = B+1)
    
    ## get original test statistics
    T_obs <- rowSums(X)/sqrt(n)
    T[, B+1] <- T_obs
    
    ## get permutation test statistic
    for (bb in 1:B) {
        eps <- rbinom(n, 1, prob = 0.5)*2 - 1  ## signs
        eX <- sweep(X, MARGIN = 2, STATS = eps, FUN = `*`)
        Tb <- rowSums(eX)/sqrt(n)
        T[, bb] <- Tb
    }
    T[, B + 1] <- T_obs
    res <- list(T = T_obs, T0 = T[, -(B+1), drop = FALSE])
    if (p.value) {
        ## get m x B matrix of pvalues under the null
        ## by sorting null test statistics as proposed by Ge et al (2003)
        ## in a permutation context
        pB <- rowRanks(-abs(T))/(B+1)
        
        res$p <- pB[, B + 1]
        res$p0 <- pB[, -(B+1), drop = FALSE]
    }
    return(res)
}


#' Perform permutation testing for a two-group model
#' 
#' @param X a matrix of \code{m} variables by \code{n} observations
#' 
#' @param cls a vector of length \code{n} class labels in \code{0,1}
#'   
#' @param B a numeric values, the number of permutations to be performed
#'   
#' @param seed an integer (or NULL) value used as a seed for random number 
#'   generation. If \code{NULL}, no seed is specified.
#'   
#' @references Ge, Y., Dudoit, S. and Speed, T.P., 2003. Resampling-based
#'   multiple testing for microarray data analysis. _Test_, 12(1), pp.1-77.
#'   
#' @return a list with two elements: \describe{
#' \item{p}{A vector of \eqn{m} \eqn{p}-values}
#' \item{p0}{A \eqn{m \times B} vector of permutation \eqn{p}-values}
#' }
#' 
#' @examples
#' 
#' p <- 513
#' n <- 38
#' mat <- matrix(rnorm(p*n), ncol=n)
#' cls <- rep(c(0, 1), times=c(27, n-27))
#' resPerm <- twoGroupPermutationTest(X=mat, cls=cls, B=1000, seed=123)
#' 
#' alpha <- 0.05
#' res <- jointFWERControl(-resPerm$p, refFamily="Simes", alpha=alpha, stat=-resPerm$p0)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
twoGroupPermutationTest <- function(X, cls, B, seed=NULL){
    ## sanity checks
    n <- ncol(X)
    if (length(cls) != n) {
        stop("The number of columns of argument 'X' should match the length of argument 'cls'")
    }
    if (!all(sort(unique(cls))==c(0,1))) {
        stop("Argument 'cls' should contain (only) 0:s and 1:s")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    m <- nrow(X)
    T <- matrix(nrow=m, ncol=B)
    
    ## Compute the observed test stat
    
    # ## fastT if we know var.equal=FALSE
    #T_obs <- genefilter::fastT(X, which(cls==0), which(cls==1), var.equal=FALSE)$z
    T_obs <- rowWelchTests(X, categ=cls)$statistic
    
    # ## rowttests if we know var.equal=TRUE
    # T_obs <- rowttests(X, factor(cls), tstatOnly=TRUE)
    
    # ## For a data-driven approach, but longer	
    # fcls <- factor(cls)
    # T_obs <- apply(X, 1, function(j) t.test(j ~ fcls)$statistic)
    
    ## get permutation test statistic
    for (bb in 1:B){
        cls_perm <- sample(cls, length(cls))
        Tb <- rowWelchTests(X, cls_perm)$statistic
        T[, bb] <- Tb
    }
    
    ## get vector of m permutation p-values
    sw <- sweep(abs(T), MARGIN=1, STATS=abs(T_obs), FUN=">=")
    p <- rowMeans(sw)
    
    ## get m x B matrix of pvalues under the null
    ## by sorting null test statistics as proposed by Ge et al (2003)
    pB <- rowRanks(-abs(T))/B
    return(list(p=p, p0=pB))
}


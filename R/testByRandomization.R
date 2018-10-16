#' Randomization-based testing
#'
#' Randomization-based testing using permutation or sign-flipping
#'
#' @param X a matrix of \code{m} variables by \code{n} observations
#'
#' @param B A numeric value, the number of permutations to be performed
#'
#' @param categ A vector of length \code{n} class labels for flavor
#'   "perm". Defaults to colnames(X).
#'   
#' @param refCat A character value, the class label to be used as a reference.
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#'
#' @param rand.p.value A boolean value: should randomization \eqn{p}-values be
#'   calculated and returned? Defaults to @FALSE
#'
#' @param seed An integer (or NULL) value used as a seed for random number
#'   generation. If \code{NULL}, no seed is specified
#'
#' @details The type of randomization is determined by the parameter \code{categ}.
#'   If \code{categ} does not contain two distinct values (or is \code{NULL}), a
#'   one-sample test is performed using randomization (flavor "flip"). If it
#'   contains two distinct values, a two-sample test is perfomed using
#'   permutations (flavor "perm").
#'
#'   For permutation, we test the null hypothesis: "both groups have the same
#'   mean" against the alternative specified by parameter \code{alternative}.
#'   The test is Welch's two-sample test for unequal variances. Permuted test
#'   statistics are calculated by B permutations of the group labels.
#'   Corresponding observed and permuted p-values are calculated as the
#'   proportion of permutations (including the identity) for which the permuted
#'   test statistic is larger than the observed test statistic.
#'
#'   For sign-flipping, we test the null hypothesis: "the mean is 0" against the
#'   alternative specified by parameter \code{alternative}. We use the
#'   (rescaled) empirical mean of the observations as a test statistic.
#'   Sign-flipped test statistics are calculated by flipping the sign of each
#'   observation with probability 1/2.
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
#'   \item{p}{A vector of \eqn{m} parametric \eqn{p}-values}
#'
#'   \item{p0}{A \eqn{m \times B} matrix of parametric \eqn{p}-values on
#'   randomized data}
#'
#'   \item{flavor}{A character value, the type of randomization performed:
#'   "perm" for permutation-based randomization in two-sample tests, and "flip"
#'   for sign-flipping-based randomization in one sample tests. See Details.}
#'
#'   \item{rand.p}{A vector of \eqn{m} \eqn{p}-values (only if
#'   \code{rand.p.value} is \code{TRUE} )}
#'
#'   \item{rand}{A \eqn{m \times B} matrix of randomization \eqn{p}-values (only
#'   if \code{rand.p.value} is \code{TRUE} )}
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
#' m <- 123
#' rho <- 0.2
#' n <- 100
#' pi0 <- 0.5
#' B <- 1e3
#'
#' ## two-sample data
#' sim <- gaussianSamples(m, rho, n, pi0, SNR = 2, prob = 0.5)
#' tests <- testByRandomization(sim$X, B, categ = colnames(sim$X))
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col = colStat, main = "Test statistics", pch  =pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#'
#' sim <- gaussianSamples(m, rho, n, pi0, SNR=2)
#' tests <- testByRandomization(sim$X, B, categ = colnames(sim$X))
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col = colStat, main = "Test statistics", pch = pch)
#' legend("topleft", c("H0", "H1"), pch = pch, col = 1:2)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
testByRandomization <- function(X, B, categ, refCat = levels(as.factor(categ))[1], 
                                alternative = c("two.sided", "less", "greater"),
                                rowTestFUN = rowWelchTests,
                                rand.p.value = FALSE, seed = NULL){
    alternative <- match.arg(alternative)
    ## sanity checks
    n <- ncol(X)
    if (missing(categ)) {
        flavor <- "flip"
        warning("No class labels ('categ') provided: performing one-sample tests by sign-flipping")
    } else {
        categ <- as.factor(categ)
        cats <- levels(categ)
        luc <- length(cats)
        if (luc <= 1) {  
            # no classes or a single class given: assuming sign flipping 
            flavor <- "flip"
        } else {
            if (length(categ) != n) {
                stop("The number of columns of argument 'X' should match the length of argument 'categ'")
            }
            if (luc == 2) {
                flavor <- "perm"
                tbl <- table(categ)
                # if ( !all(names(tbl) == c("0", "1"))) {  # note that numeric values are allowed in categ as they are converted into character by 'table'...
                #     stop("Argument 'categ' should contain (only) 0:s and 1:s")
                # }
                if (min(tbl) < 3) {
                    stop("Argument 'categ' should contain at least 3 elements of each sample")
                }
            } else if (luc > 2) {
                stop("Tests for more than 2 classes not implemented yet")
            }
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
        
        ## observed
##        rwt <- rowWelchTests(X, categ = categ, refCat = refCat, alternative = alternative)
        rwt <- rowTestFUN(X, categ = categ, refCat = refCat, alternative = alternative)
        T <- rwt$statistic
        p <- rwt$p.value  ## parametric p-value
        df <- rwt$parameter  ## degrees of freedom of the T statistics
        rm(rwt)
        
        ## under H0
        T0 <- matrix(nrow = m, ncol = B) ## test statistics under the null
        p0 <- matrix(nrow = m, ncol = B) ## parametric p-value
        df0 <- matrix(nrow = m, ncol = B) 
        for (bb in 1:B) {
            categ_perm <- sample(categ, length(categ))
            ## rwt <- rowWelchTests(X, categ = categ_perm, refCat = refCat, alternative = alternative)
            rwt <- rowTestFUN(X, categ = categ_perm, refCat = refCat, alternative = alternative)
            T0[, bb] <- rwt$statistic
            p0[, bb] <- rwt$p.value
            df0[, bb] <- rwt$parameter
        }
        res <- list(T = T, T0 = T0, 
                    flavor = flavor,
                    p = p, p0 = p0,
                    df = df, df0 = df0)
    } else if (flavor == "flip") {
        ## observed test statistics and p-values
        T <- rowSums(X)/sqrt(n)
        p <- switch(alternative, 
                        "two.sided" = 2*(1 - pnorm(abs(T))),
                        "greater" = 1 - pnorm(T),
                        "less" = pnorm(T))
        ## test statistics under H0
        T0 <- testBySignFlipping(X, B)
        p0 <- switch(alternative, 
                     "two.sided" = 2*(1 - pnorm(abs(T0))),
                     "greater" = 1 - pnorm(T0),
                     "less" = pnorm(T0))
        res <- list(T = T, T0 = T0, p = p, p0 = p0, flavor = flavor)
    }
    
    if (rand.p.value) {
        ## get m x (B+1) matrix of pvalues under the null (+ original)
        ## by sorting null test statistics as proposed by Ge et al (2003)
        TT <- cbind(T0, T)
        pB <- switch(alternative, 
                       "two.sided" = rowRanks(-abs(TT)) / (B+1),
                       "greater" = rowRanks(-TT) / (B+1),
                       "less" = rowRanks(TT) / (B+1))

        res$rand.p <- pB[, B+1]
        res$rand.p0 <- pB[, -(B+1), drop = FALSE]
    }
    return(res)
}

# not used!
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

# not used!
testByPermutationR <- function(X, categ, B) {
    m <- nrow(X)
    n <- ncol(X)
    stopifnot(n == length(categ))
    
    T <- matrix(nrow = m, ncol = B)
    for (bb in 1:B) {
        categ_perm <- sample(categ, length(categ))
        Tb <- rowWelchTests(X, categ = categ_perm)$statistic
        T[, bb] <- Tb
    }
    T
}
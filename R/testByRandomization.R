#' Randomization-based testing
#'
#' Randomization-based testing using permutation or sign-flipping
#'
#' @param X a matrix of \code{m} variables by \code{n} observations
#'
#' @param categ An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#'   
#' @param B A numeric value, the number of permutations to be performed
#'
#' @param rowTestFUN A (vectorized) test function used in the two-sample case.
#'   Defaults to \code{\link{rowWelchTests}}
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#'
#' @param rand.p.value A boolean value: should randomization \eqn{p}-values be
#'   calculated and returned? Defaults to @FALSE
#'
#' @details The type of randomization is determined by the column names of
#'   \code{X}. If these column names have exactly two distinct values, the
#'   corresponding columns are interpreted as two samples and a two-sample
#'   permutation-based test  is performed (flavor "perm"). Otherwise (including
#'   if \code{X} does not have column names), a one-sample test is performed
#'   using sign-flipping (flavor "flip").
#'
#'   For permutation, we test the null hypothesis: "both groups have the same
#'   mean" against the alternative specified by parameter \code{alternative}. By
#'   default, the test is Welch's two-sample test for unequal variances, but
#'   other tests may be used via the argument \code{rowTestFUN}.  Permuted test
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
#'}
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
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
#' X <- sim$X
#' categ <-sim$categ
#' tests <- testByRandomization(X = X, categ = categ, B)
#' testsW <- testByRandomization(X = X, categ = categ, B = 10, rowTestFUN = rowWilcoxonTests)
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col = colStat, main = "T-Test statistics", pch = pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#'
#' plot(testsW$T, col = colStat, main = "Wilcoxon test statistics", pch = pch)
#' legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
#' 
#' # one-sample data:
#' sim <- gaussianSamples(m, rho, n, pi0, SNR=2)
#' tests <- testByRandomization(X = sim$X, categ = sim$categ, B, alternative = "two.sided")
#'
#' ## show test statistics
#' pch <- 20
#' colStat <- 1+sim$H
#' plot(tests$T, col = colStat, main = "Test statistics", pch = pch)
#' legend("topleft", c("H0", "H1"), pch = pch, col = 1:2)
#'
#' plot(-log10(tests$p), col = colStat, main = "-log[10](p-value)", pch = pch)
#' legend("topleft", c("H0", "H1"), pch = pch, col = 1:2)
#'
#' @importFrom matrixStats rowRanks
#' @export
#' 
testByRandomization <- function(X, categ, B, 
                                alternative = c("two.sided", "less", "greater"),
                                rowTestFUN = rowWelchTests,
                                rand.p.value = FALSE){
    alternative <- match.arg(alternative)
    ## sanity checks
    n <- ncol(X)
    if (missing(categ) || (length(unique(categ))==1)) {
        flavor <- "flip"
    } else {
        if (length(unique(categ)) == 2) {
            stopifnot(all(categ %in% c(0,1)))
            flavor <- "perm"
            tbl <- table(categ)
            if (min(tbl) < 3) {
                stop("At least 3 observations from each sample are required for two-sample tests")
            }
        } else {
            stop("Tests for more than 2 classes not implemented yet")
        }
    }
    m <- nrow(X)
    
    if (flavor == "perm") {
        ## TODO: (cf issue #3)
        ## * other statistics ? (difference in empirical means, Mann-Whitney)
        ## * one-sided tests ?
        
        ## map class labels to 0-1 for simplicity of implementation of rowTestFUN(s)
        levels(categ) <- c("0", "1") 
        ## observed
        rwt <- rowTestFUN(X, categ = categ, alternative = alternative)
        T <- rwt$statistic
        p <- rwt$p.value  ## parametric p-value
        # df <- rwt$parameter  ## degrees of freedom (for T tests; possibly NULL for other tests)
        rm(rwt)
        
        ## under H0
        T0 <- matrix(nrow = m, ncol = B) ## test statistics under the null
        p0 <- matrix(nrow = m, ncol = B) ## parametric p-value
        # df0 <- matrix(nrow = m, ncol = B) 
        for (bb in 1:B) {
            categ_perm <- sample(categ, length(categ))
            rwt <- rowTestFUN(X, categ = categ_perm, alternative = alternative)
            T0[, bb] <- rwt$statistic
            p0[, bb] <- rwt$p.value
        }
        res <- list(T = T, T0 = T0, 
                    flavor = flavor,
                    p = p, p0 = p0)
                    # df = df, df0 = df0)
    } else if (flavor == "flip") {
        ## observed test statistics and p-values
        T <- rowSums(X)/sqrt(n)
        p <- switch(alternative, 
                    "two.sided" = 2*(pnorm(abs(T), lower.tail = FALSE)),
                    "greater" = pnorm(T, lower.tail = FALSE),
                    "less" = pnorm(T, lower.tail = TRUE))
        ## test statistics under H0
        T0 <- testBySignFlipping(X, B)
        p0 <- switch(alternative, 
                     "two.sided" = 2*(pnorm(abs(T0), lower.tail = FALSE)),
                     "greater" = pnorm(T0, lower.tail = FALSE),
                     "less" = pnorm(T0, lower.tail = TRUE))
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
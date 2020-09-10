# Empirical coverage of a threshold family
# 
# Empirical coverage of a threshold family
#
# @param thr A numeric vector of length \code{m}, a threshold family. Should  be sorted in increasing order.
# @param mat A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
# samples of \code{c} p-values *under the null hypothesis*.
# 
# @examples
# alpha <- 0.05
# 
# # empirical coverage of Simes thresholds for independent Gaussian nulls
# m <- 234
# B <- 1234
# X0 <- replicate(B, rnorm(m))
# p0 <- pnorm(X0, lower.tail = FALSE)
# thr <- SimesThresholdFamily(m)(alpha)
# sansSouci:::empiricalCoverage(thr, p0)
# 
# @return The empirical coverage of \code{thr} with respect to \code{mat}.
empiricalCoverage <- function(thr, mat) {
    kMax <- length(thr)
    stopifnot(kMax<=nrow(mat)) ## sanity check
    kmaxH0 <- -partialColSortDesc(-mat, kMax); ## Implicitly forcing 'Rcpp' flavor
    empiricalCoverageO(-thr, -kmaxH0);
}


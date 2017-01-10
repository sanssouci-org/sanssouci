#' Empirical coverage of a threshold family
#' 
#' Empirical coverage of a threshold family
#'
#' @param thr A numeric vector of length \code{m}, a threshold family. Should  be sorted in decreasing order.
#' @param mat A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
#' samples of \code{c} test statistics under the null hypothesis.

#' @return The empirical coverage of \code{thr} with respect to \code{mat}.
empiricalCoverage <- function(thr, mat) {
    kMax <- length(thr)
    stopifnot(kMax<=nrow(mat)) ## sanity check
    kmaxH0 <- partialColSortDesc(mat, kMax); ## Implicitly forcing 'Rcpp' flavor
    empiricalCoverageO(thr, kmaxH0);
}


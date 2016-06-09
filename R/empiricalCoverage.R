##' @export
empiricalCoverage <- function(thr, mat) {
    kMax <- length(thr)
    stopifnot(kMax<=nrow(mat)) ## sanity check
    kmaxH0 <- partialColSortDesc(mat, kMax); ## Implicitly forcing 'Rcpp' flavor
    empiricalCoverageO(thr, kmaxH0);
}


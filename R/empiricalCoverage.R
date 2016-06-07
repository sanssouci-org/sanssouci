##' @export
empiricalCoverage <- function(thr, mat, kMax=nrow(mat)) {
    kmaxH0 <- partialColSortDesc(mat, kMax); ## Implicitly forcing 'Rcpp' flavor
    empiricalCoverageO(thr, kmaxH0);
}


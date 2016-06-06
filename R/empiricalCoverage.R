empiricalCoverage <- function(thr, mat, kMax=ncol(mat)) {
    kmaxH0 <- partialColSortDesc(mat, kMax); ## Implicitly forcing 'Rcpp' flavor
    empiricalCoverageO(thr, kmaxH0);
}


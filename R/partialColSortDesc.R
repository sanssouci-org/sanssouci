# partialColSortDesc
#
# partial sorting of the columns of a matrix by descending order
#
#
# @param mat a \code{m} x \code{B} numeric matrix
# @param k an integer value between 1 and \code{m}
# @param Rcpp sorting is performed in C++ if \code{TRUE} (the default), and
# using \code{base::sort} if \code{FALSE}.
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
# @return a \code{k} x \code{B} numeric matrix, whose \code{b}-th
# column corresponds to the \code{k} largest values of the
# \code{b}-th column of \code{mat}, sorted decreasingly.
# @examples
#
# A <- matrix(rnorm(15), 5, 3);
# print(A)
# BR <- partialColSortDesc(A, 2, Rcpp=FALSE);
# B <- partialColSortDesc(A, 2, Rcpp=TRUE);
# stopifnot(identical(B,BR));
#
# if (require("rbenchmark")) {
#     A <- matrix(rnorm(1e6), 10000, 1000);
#     kmax <- 100;
#     benchmark(partialColSortDesc(A, kmax, Rcpp=TRUE),
#               partialColSortDesc(A, kmax, Rcpp=FALSE),
#               replications=1)
# }
#
partialColSortDesc <- function(
    mat,
    k=nrow(mat),
    Rcpp=TRUE)   {
    k <- round(k)  ## make sure 'k' is not converted to the wrong integer as in 'as.integer(1000*(1-0.8))==199'
    stopifnot( 1 <= k && k <= nrow(mat))
    if (Rcpp) {
        pcsd <- partialColSortDescCpp(mat, k);
    } else {
        pcsd <- -apply(-mat, 2, sort, partial=1:k)
        pcsd <- pcsd[1:k, , drop=FALSE]
    }
    ## sanity checks
    stopifnot(nrow(pcsd)==k)
    stopifnot(ncol(pcsd)==ncol(mat))
    return(pcsd)
}

partialColSort <- function(
    mat,
    k=nrow(mat),
    Rcpp=TRUE)   {
    k <- round(k)  ## make sure 'k' is not converted to the wrong integer as in 'as.integer(1000*(1-0.8))==199'
    stopifnot( 1 <= k && k <= nrow(mat))
    if (Rcpp) {
        pcsd <- partialColSortCpp(mat, k);
    } else {
        pcsd <- apply(mat, 2, sort, partial=1:k)
        pcsd <- pcsd[1:k, , drop=FALSE]
    }
    ## sanity checks
    stopifnot(nrow(pcsd)==k)
    stopifnot(ncol(pcsd)==ncol(mat))
    return(pcsd)
}

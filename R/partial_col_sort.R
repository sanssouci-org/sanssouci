#' partial_col_sort
#'
#' partial sorting of the columns of a matrix
#'
#' @param mat a \code{m} x \code{B} numeric matrix
#' @param k an integer value between 1 and \code{m}
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#'
#' @return a \codée{k} x \code{B} numeric matrix, whose \code{b}-th
#' column corresponds to the \code{k} largest values of the
#' \code{b}-th column of \code{mat}, sorted decreasingly.
#'
#' @noRd
#'
#' @examples
#'
#' A <- matrix(rnorm(15), 5, 3)
#' print(A)
#' partial_col_sort(A, 2)
#'
partial_col_sort <- function(mat, k = nrow(mat)) {
  k <- round(k) ## make sure 'k' is not converted to the wrong integer
  # as in 'as.integer(1000*(1-0.8))==199'
  stopifnot(1 <= k && k <= nrow(mat))
  pcsd <- apply(mat, 2, sort, partial = 1:k)
  pcsd <- pcsd[1:k, , drop = FALSE]

  ## sanity checks
  stopifnot(nrow(pcsd) == k)
  stopifnot(ncol(pcsd) == ncol(mat))
  return(pcsd)
}

#' @inherit partial_col_sort
#' @noRd
#' @examples
#' A <- matrix(rnorm(15), 5, 3)
#' partial_col_sort_desc(A, 2)
partial_col_sort_desc <- function(mat, k = nrow(mat)) {
  k <- round(k) ## make sure 'k' is not converted to the wrong
  ## integer as in 'as.integer(1000*(1-0.8))==199'
  stopifnot(1 <= k && k <= nrow(mat))
  pcsd <- -apply(-mat, 2, sort, partial = 1:k)
  pcsd <- pcsd[1:k, , drop = FALSE]

  ## sanity checks
  stopifnot(nrow(pcsd) == k)
  stopifnot(ncol(pcsd) == ncol(mat))
  return(pcsd)
}

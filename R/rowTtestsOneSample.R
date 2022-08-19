#' One sample T tests for rows of a matrix
#'
#' @param mat A numeric matrix whose rows correspond to \code{m} variables and columns to \code{n}
#'   observations
#'
#' @param labels Either a numeric vector of \code{n} values in \eqn{-1, 1} to
#'   perform sign flips, or a \code{n x B} matrix stacking \code{B} such
#'   vectors. If missing, coerced to \code{rep(1, n)}
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". `Alternative = 
#'   "greater"` corresponds to a positive mean.
#'
#' @return A list containing the following components:
#' \describe{ 
#'   \item{statistic}{the value of the statistics}
#'   \item{p.value}{the p-values for the tests}}
#'   Each of these elements is a matrix of size \code{nrow(mat) x B}, coerced to a vector of length \code{nrow(mat)} if \code{B=1}
#' @author Pierre Neuvial & Nicolas Enjalbert Courrech 
#'

rowTtestsOneSample <- function(mat, labels, alternative = c("two.sided", "grater", "less"), mu=0) {
  n <- ncol(mat)
  if (missing(labels)) {
    labels <- rep(1, n)
  } 
  if(is.null(dim(labels))){
    B <- 1
  } else {
    B <- dim(labels)[2]
  }
  mean_ <- mat %*% labels / n
  var_ <- ((rowSums(mat^2) %*% t(rep(1,B))) - n * mean_^2 )/(n-1)
  test <- ((mean_ - mu )/ sqrt(var_) ) * sqrt(n)
  df <- n-1
  pval <- switch(alternative,
                 "two.sided" = 2*(1 - pt(abs(test), df = df)),
                 "greater" = 1 - pt(test, df = df),
                 "less" = pt(test, df = df))
  list(statistic = test, p.value = pval)
}

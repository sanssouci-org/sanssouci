#' Pearson's correlation test for rows of a matrix
#'
#' @param X A \code{m x n} numeric matrix whose rows correspond to variables
#'   and columns to observations
#'
#' @param groups Either a numeric vector of continuous covariate for
#'   the observations
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{t.test}}, alternative = "greater" is the alternative that class
#'   1 has a larger mean than class 0.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{statistic}{the value of the t-statistics}
#'   \item{p_value}{the p-values for the tests}
#'   \item{estimate}{the correlation}
#'   \item{parameter}{the degrees of freedom for the t-statistics}}
#'   Each of these elements is a matrix of size \code{m x B},
#'   coerced to a vector of length \code{m} if \code{B=1}
#'
#' @details This function is a wrapper around the
#' \code{row_cor_pearson function} in the 'matrixTests' package.
#'
#' @seealso [matrixTests::row_cor_pearson()]
#'
#' @importFrom matrixTests row_cor_pearson
#'
#' @author Nicolas Enjalbert Courrech
#'
#' @export
#'
#' @examples
#'
#' m <- 300
#' n <- 38
#' mat <- matrix(rnorm(m * n), ncol = n)
#' groups <- rnorm(n, mean = 10)
#' system.time(fwt <- row_pearson_correlation_tests(mat, groups,
#'   alternative = "greater"
#' ))
#' str(fwt)
#'
row_pearson_correlation_tests <- function(X, groups, alternative = "two.sided") {
  if (is.null(dim(groups))) {
    r <- row_cor_pearson(X, groups, alternative = alternative)
    r$p_value <- r$pvalue
  } else {
    m <- dim(X)[1]
    B <- dim(groups)[2]
    r <- list(
      p_value = matrix(nrow = m, ncol = B),
      statistic = matrix(nrow = m, ncol = B),
      cor = matrix(nrow = m, ncol = B)
    )
    for (b in 1:B) {
      rr <- row_cor_pearson(X, groups[, b], alternative = alternative)
      r$p_value[, b] <- rr$pvalue
      r$statistic[, b] <- rr$statistic
      r$cor[, b] <- rr$cor
    }
  }
  return(list(
    statistic = r$statistic,
    p_value = r$p_value,
    estimate = r$cor
  ))
}

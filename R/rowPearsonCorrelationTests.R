#' Pearson's correlation test for rows of a matrix
#'
#' @param X A \code{m x n} numeric matrix whose rows correspond to variables
#'   and columns to observations
#'
#' @param categ Either a numeric vector of continuous covariate for
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
#'   \item{parameter}{the degrees of freedom for the t-statistics}
#'   \item{p.value}{the p-values for the tests} 
#'   \item{estimate}{the correlation}}
#'   Each of these elements is a matrix of size \code{m x B}, coerced to a vector of length \code{m} if \code{B=1}
#'   
#' @details This function is a wrapper around the \code{row_cor_pearson function} in the 'matrixTests' package.
#' 
#' @seealso [matrixTests::function(row_cor_pearson)]
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
#' mat <- matrix(rnorm(m*n), ncol=n)
#' categ <- rnorm(n, mean = 10)
#' system.time(fwt <- rowPearsonCorrelationTests(mat, categ, alternative = "greater"))
#' str(fwt)
#' 
rowPearsonCorrelationTests <- function(X,categ,alternative = "two.sided"){
  if(is.null(dim(categ))){
    r <- row_cor_pearson(X, categ, alternative = alternative)
  } else {
    m <- dim(X)[1]
    B <- dim(categ)[2]
    r <- list(pvalue = matrix(nrow = m, ncol = B), statistic = matrix(nrow = m, ncol= B), cor = matrix(nrow = m, ncol = B))
    for (b in 1:B){
      rr <- row_cor_pearson(X, categ[,b], alternative = alternative)
      r$pvalue[,b] <- rr$pvalue
      r$statistic[,b] <- rr$statistic
      r$cor[,b] <- rr$cor
    }
  }
  return(list(p.value = r$pvalue, statistic = r$statistic, estimate =r$cor))
}
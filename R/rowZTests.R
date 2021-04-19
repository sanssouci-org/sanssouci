#' Z tests for rows of a matrix
#'
#' @param mat A numeric matrix whose rows correspond to \code{m} variables and columns to \code{n}
#'   observations
#'
#' @param categ Either a numeric vector of \code{n} values in \eqn{-1, 1} to
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
#' @author Pierre Neuvial
#'
#' @export
#'
#' @examples
#'
#' p <- 1e4+1
#' n <- 380
#' mat <- matrix(rnorm(p*n, mean = 1), ncol=n)
#' zt <- rowZTests(mat, alternative = "greater")
#' str(zt)
#'
#' # compare with apply version:
#' p <- apply(mat, 1, FUN=function(x) {
#'    stat <- sum(x)/sqrt(length(x))
#'    pnorm(stat, lower.tail = FALSE)
#' })
#' all(abs(zt$p.value - p) < 1e-10)  ## same results
#' 
#' # Sign flipping
#' B <- 10
#' eps <- replicate(B, rbinom(n, 1, 0.5)*2 - 1)  ## Rademacher
#' zt_perm <- rowZTests(mat, eps, alternative = "greater")
#' str(zt_perm)

rowZTests <- function(mat, categ, alternative) {
    n <- ncol(mat)
    if (missing(categ)) {
        categ <- rep(1, n)
    } 
    T <- mat %*% categ / sqrt(n)
    p <- switch(alternative, 
                "two.sided" = 2*(pnorm(abs(T), lower.tail = FALSE)),
                "greater" = pnorm(T, lower.tail = FALSE),
                "less" = pnorm(T, lower.tail = TRUE))
    list(statistic = T, p.value = p)
}

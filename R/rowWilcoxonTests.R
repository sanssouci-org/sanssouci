#' Wilcoxon rank sum tests for each row of a matrix
#'
#' Vectorized version of two-sample Wilcoxon rank sum tests
#'
#' @param mat A numeric matrix whose rows correspond to variables and columns to
#'   observations
#'
#' @param categ A numeric vector of \code{ncol(mat)} categories in \eqn{0, 1} for the
#'   observations
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{wilcox.test}}, alternative = "greater" is the alternative that
#'   class 1 is shifted to the right of class 0.
#'
#' @param correct A logical indicating whether to apply continuity correction in
#'   the normal approximation for the p-value.
#'
#' @return A list with class "htest" containing the following components:
#'   \describe{ \item{statistic}{the value of the statistics} \item{p.value}{the
#'   p-values for the tests}} 
#'
#' @details The p-values are computed using the normal approximation as
#'   described in the \code{\link{wilcox.test}} function. The exact p-values
#'   (which can be useful for small samples with no ties) are not implemented
#'   (yet).
#'
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @seealso wilcox.test
#' @return A data.frame with columns \describe{ \item{stat}{A vector of \code{m}
#'   Wilcoxon sum rank test statistics of association between \code{X} and
#'   \code{y}.} \item{stat0Mat}{An \code{m} x \code{B} matrix of \code{B}
#'   realizations of a \code{m}-dimensional vector of test statistics under the
#'   null hypothesis of no association between \code{X} and \code{y}.}}
#' @importFrom matrixStats rowRanks rowTabulates
#' @export
#' @examples
#'
#' p <- 100
#' n <- 120
#' mat <- matrix(rnorm(p*n), ncol = n)
#' cls <- rep(c(0, 1), times = c(n/2, n/2))
#' fwt <- rowWilcoxonTests(mat, categ = cls, alternative = "two.sided")
#' str(fwt)
#' 
#' # compare with ordinary wilcox.test:
#' pwt <- t(sapply(1:p, FUN=function(ii) {
#'   wt <- wilcox.test(mat[ii, cls==1], mat[ii, cls==0], alternative = "two.sided")
#'   c(statistic = wt[["statistic"]], p.value = wt[["p.value"]])
#' }))
#' all(abs(fwt$p.value-pwt[, "p.value"]) < 1e-10)  ## same results
#' all(abs(fwt$statistic-pwt[, "statistic.W"]) < 1e-10)  ## same results
#' 
rowWilcoxonTests <- function(mat, categ, alternative = c("two.sided", "less", "greater"), correct = TRUE) {
    alternative <- match.arg(alternative)
    categCheck(categ, ncol(mat))

    rks <- rowRanks(mat, ties.method = "average")
    ny <- sum(categ == 0)
    nx <- sum(categ == 1)
    min_stat <- nx*(nx + 1)/2  ## mininmal rank sum

    wx <- which(categ == 1)
    stat <- rowSums(rks[, wx])
    stat <- stat - min_stat
    
    ## gaussian approximation for p-values (presence of ties or n > 50)
    ## source: 'stats:::wilcox.test.default'
    mode(rks) <- "integer"   # rowTabulates only takes integer values
    n_ties <- rowTabulates(rks)  
    ties <- rowSums(n_ties^3 - n_ties) 
    sigma <- sqrt((nx * ny / 12) * ((ny + nx + 1) - ties/((ny + nx) * (ny + nx - 1))))
    
    z <- stat - ny*nx/2
    CORRECTION <- 0
    if (correct) {
        CORRECTION <- switch(alternative, two.sided = sign(z) * 
                                 0.5, greater = 0.5, less = -0.5)
    }
    z <- (z - CORRECTION)/sigma
    p <- switch(alternative, 
                less = pnorm(z), 
                greater = pnorm(z, lower.tail = FALSE), 
                two.sided = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE)))

    data.frame(statistic = stat,
                p.value = p)
}

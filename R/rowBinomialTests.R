#' Binomial proportion tests for each row of a matrix
#'
#' Vectorized version of the binomial test
#'
#' @param mat A numeric matrix whose rows correspond to variables and columns to
#'   observations
#'
#' @param categ A numeric vector of \code{ncol(mat)} categories in \eqn{0,1} for the
#'   observations
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{binom.test}}, alternative = "greater" is the alternative that
#'   class 1 is shifted to the right of class 0.
#'
#' @param warn A boolean value indicating whether to issue a warning if
#'   \code{alternative=="two-sided"}. Defaults to \code{TRUE}.
#'
#' @return A list with class "htest" containing the following components:
#'   \describe{ \item{statistic}{the value of the statistics} \item{p.value}{the
#'   p-values for the tests}}
#'
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @seealso binom.test
#' @return A data.frame with columns \describe{ \item{stat}{A vector of \code{m}
#'   Wilcoxon sum rank test statistics of association between \code{X} and
#'   \code{y}.} \item{stat0Mat}{An \code{m} x \code{B} matrix of \code{B}
#'   realizations of a \code{m}-dimensional vector of test statistics under the
#'   null hypothesis of no association between \code{X} and \code{y}.}}
#' @importFrom matrixStats rowRanks rowTabulates
#' @export
#' @examples
#'
#' alt <- c("two.sided", "less", "greater")[1]
#'
#' p <- 100
#' n0 <- 60; n1 <- 40
#' mat0 <- matrix(rbinom(p*n0, size = 1, prob = 0.05), ncol = n0)
#' mat1 <- matrix(rbinom(p*n1, size = 1, prob = 0.02), ncol = n1)
#' mat <- cbind(mat0, mat1)
#' cls <- rep(c(0, 1), times = c(n0, n1))
#' fbt <- rowBinomialTests(mat, categ = cls, alternative = alt)
#' str(fbt)
#'
#' # compare with ordinary binom.test:
#' pbt <- t(sapply(1:p, FUN=function(ii) {
#'   x1 <- mat[ii, cls==1]
#'   x0 <- mat[ii, cls==0]
#'   bt <- binom.test(sum(x1), length(x1), mean(x0), alternative = alt)
#'   c(statistic = bt[["statistic"]], p.value = bt[["p.value"]])
#' }))
#' all(abs(fbt$p.value-pbt[, "p.value"]) < 1e-10)  ## same results
#' all(abs(fbt$statistic-pbt[, "statistic.number of successes"]) < 1e-10)  ## same results
#' 
#' @importFrom stats pbinom dbinom
rowBinomialTests <- function(mat, categ, alternative = c("two.sided", "less", "greater"), warn = TRUE) {
    alternative <- match.arg(alternative)
    if (alternative == "two.sided" & warn) {
        warning("Two-sided p-value not vectorized yet! Looping for now.")
    }
    categCheck(categ, ncol(mat))

    nx <- sum(categ == 1)
    x <- rowSums(mat[, categ == 1])
    py <- rowMeans(mat[, categ == 0])  # parameter of the control distribution
    p <- switch(alternative,
                less = pbinom(x, nx, py), 
                greater = pbinom(x - 1, nx, py, lower.tail = FALSE), 
                two.sided = rowBinomialP_2s(x = x, n = nx, p = py, relErr = 1 + 1e-07))
    data.frame(statistic = x,
                p.value = p)
}

# taken from the source code of stats::binomial.test
rowBinomialP_2s <- function(x, n, p, relErr = 1 + 1e-07) {
    #stopifnot(length(x) == length(p) || (length(x) == 1) || length(p) == 1)
    stopifnot(length(x) == length(p))
    pval <- rep(NA_real_, length(x))
    pval[which(p == 0)] <- (x[which(p == 0)] == 0)
    pval[which(p == 1)] <- (x[which(p == 1)] == n)
    ww <- which(p > 0 | p < 1)
    for (jj in ww) { 
        xx <- x[jj]
        pp <- p[jj]
        d <- dbinom(xx, n, pp)
        m <- n * pp
        if (xx == m) {
            pval[jj] <- 1
        } else if (xx < m) {
            i <- seq.int(from = ceiling(m), to = n)
            y <- sum(dbinom(i, n, pp) <= d * relErr)
            pval[jj] <- pbinom(xx, n, pp) + pbinom(n - y, n, pp, lower.tail = FALSE)
        } else {
            i <- seq.int(from = 0, to = floor(m))
            y <- sum(dbinom(i, n, pp) <= d * relErr)
            pval[jj] <- pbinom(y - 1, n, pp) + pbinom(xx - 1, n, pp, lower.tail = FALSE)
        }
    }
    pval
}
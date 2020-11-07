#' Welch t-tests for rows of a matrix
#'
#' Welch t-tests for each row of a matrix, intended to be speed efficient
#'
#' Other functions to perform t-tests on each row of a matrix include:
#' \code{\link[genefilter]{fastT}}: this function only calculates the test
#' statistics and not the corresponding degrees of freedom or p-value;
#' \code{\link[genefilter]{rowttests}}: this function calculates the test
#' statistics, degrees of freedom and associated p-value, but only in the case
#' of equal group variances (standard t-test and not Welch t-test)
#'
#' @param mat A numeric matrix whose rows correspond to variables and columns to
#'   observations
#'
#' @param categ A numeric vector of \code{ncol(mat)} categories in \eqn{0, 1} for the
#'   observations
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{t.test}}, alternative = "greater" is the alternative that
#'   class 1 has a larger mean than class 0.
#'
#' @return A data.frame with containing the following columns:
#'   \describe{ \item{statistic}{the value of the t-statistics}
#'   \item{parameter}{the degrees of freedom for the t-statistics}
#'   \item{p.value}{the p-values for the tests}} \item{meanDiff}{the mean
#'   difference}
#' @author Pierre Neuvial
#'
#' @references B. L. Welch (1951), On the comparison of several mean values: an
#'   alternative approach. Biometrika, *38*, 330-336
#'
#' @export
#'
#' @examples
#'
#' p <- 1e3
#' n <- 38
#' mat <- matrix(rnorm(p*n), ncol=n)
#' categ <- rep(c(0, 1), times=c(27, n-27))
#' fwt <- rowWelchTests(mat, categ, alternative = "greater")
#' str(fwt)
#'
#' # compare with ordinary t.test:
#' pwt <- apply(mat, 1, FUN=function(x) {
#'    t.test(x[categ==1], x[categ==0], alternative = "greater")$p.value
#' })
#' all(abs(fwt$p.value-pwt) < 1e-10)  ## same results
#' 
rowWelchTests <- function(mat, categ, alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    categCheck(categ, ncol(mat))

    sstats <- getSummaryStats(mat, categ = categ)
    Y <- sstats[["0"]]
    X <- sstats[["1"]]  ## as per the doc of t.test:
    ## 'alternative = "greater"' is the alternative that 'x' has a larger  mean
    ## than 'y'.
    swt <- suffWelchTest(X[["mean"]], Y[["mean"]],
                         X[["sd"]], Y[["sd"]],
                         X[["n"]], Y[["n"]],
                         alternative = alternative)
    swt$meanDiff <- X[["mean"]] - Y[["mean"]]
    swt
}



categCheck <- function(categ, n) {
    stopifnot(length(categ) == n)
    categ <- as.factor(categ)
    cats <- levels(categ)
    if (!identical(cats, c("0", "1"))) {
        stop("Expected two categories named '0' and '1'!")
    }
}


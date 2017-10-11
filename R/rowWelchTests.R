#' Welch t-tests for rows of a matrix
#' 
#' Welch t-tests for each row of a matrix, intended to be speed efficient
#' 
#' Other functions to perform t-tests on each row of a matrix include:
#' \code{\link[genefilter]{fastT}}: this function only calculates the test
#' statistics and not the corresponding of freedom or p-value; 
#' \code{\link[genefilter]{rowttests}}: this function calculates the test
#' statistics, degrees of freedom and associated p-value, but only in the case
#' of equal group variances (standard t-test and not Welch t-test)
#' 
#' @param mat A numeric matrix whose rows correspond to variables and columns to
#'   observations
#'   
#' @param categ A vector of \code{ncol(mat)} categories for the observations
#'   
#' @return A list with class "htest" containing the following components: 
#'   \describe{ \item{statistic}{the value of the t-statistics} 
#'   \item{parameter}{the degrees of freedom for the t-statistics} 
#'   \item{p.value}{the p-values for the tests}}
#'   
#' @author Pierre Neuvial
#'   
#' @references B. L. Welch (1951), On the comparison of several mean values: an
#'   alternative approach. Biometrika, *38*, 330-336
#'   
#' @examples
#' 
#' if (require("multtest")) { 
#'   data(golub, package="multtest")
#'   fwt <- fastWelchTest(golub, categ=golub.cl)
#'   str(fwt)
#' } else {
#'   p <- 1e5
#'   n <- 38
#'   mat <- matrix(rnorm(p*n), ncol=n)
#'   cls <- rep(c(0, 1), times=c(27, n-27))
#'   fwt <- rowWelchTests(mat, categ=cls)
#'   str(fwt)
#'   
#'   # compare with ordinary t.test:
#'   pwt <- apply(mat, 1, FUN=function(x) {
#'      t.test(x[cls==0], x[cls==1])$p.value
#'   })
#'   sum(abs(fwt$p.value-pwt))  ## same results
#' }
#' 
rowWelchTests <- function(mat, categ=colnames(mat)) {
    cats <- unique(categ)
    if (length(cats) != 2) {
        stop("Two categories expected!")
    }
    sstats <- getSummaryStats(mat, categ = categ)
    X <- sstats[[1]]
    Y <- sstats[[2]]
    swt <- suffWelchTest(X[["mean"]], Y[["mean"]],
                           X[["sd"]], Y[["sd"]],
                           X[["n"]], Y[["n"]])
    swt
}
        
#' Welch test from sufficient statistics
#'
#' @param mx A numeric value or vector, the sample average for condition "x"
#' @param my A numeric value or vector of the same length as 'mx', the sample
#'   average for condition "y"
#' @param sx A numeric value or vector of the same length as 'mx', the standard
#'   deviation for condition "x"
#' @param sy A numeric value or vector of the same length as 'mx', the standard
#'   deviation for condition "y"
#' @param nx A numeric value or vector of the same length as 'mx', the sample
#'   size for condition "x"
#' @param ny A numeric value or vector of the same length as 'mx', the sample
#'   size for condition "y"
#' @param alternative A character string specifying the alternative hypothesis.
#'   Currently only "two.sided" is implemented
#' @return A list with elements:
#'   \describe{ \item{statistic}{the value of the t-statistic}
#'   \item{parameter}{the degrees of freedom for the t-statistic}
#'   \item{p.value}{the p-value for the test} }
#' @author Pierre Neuvial
#'
#' @details In accordance with the implementation of \code{\link[stats]{t.test}} and friends, we
#'   follow the rule that 'alternative = "greater"' is the alternative that 'x'
#'   has a larger  mean than 'y'.
#'
#' @export
#' @importFrom stats pt
#' @examples
#'
#' # Reproducing the results of the 't.test' function
#'
#' x <- rnorm(100)
#' y <- rnorm(34, mean = 1)
#' target <- t.test(x, y)
#' swt <- suffWelchTest(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
#' print(swt$statistic - target$statistic)
#' print(swt$p.value - target$p.value)
#' print(swt$parameter-target$parameter)
#' 
#' target <- t.test(x, y, alternative = "greater")
#' swt <- suffWelchTest(mean(x), mean(y), sd(x), sd(y), length(x), length(y), alternative = "greater")
#' print(swt$statistic - target$statistic)
#' print(swt$p.value - target$p.value)
#' print(swt$parameter-target$parameter)
suffWelchTest <- function(mx, my, sx, sy, nx, ny, 
                            alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
#    if (alternative != "two.sided") stop("altenative", alternative, "not implemented yet!")
    
    ## sanity checks
    p <- length(mx)
    stopifnot(length(my)==p)
    stopifnot(length(sx)==p)
    stopifnot(length(sy)==p)
    stopifnot(length(nx) %in% c(1, p))
    stopifnot(length(ny) %in% c(1, p))

    sse.x <- sx^2/nx  ## sc <- sqrt((sum2c - sumc^2/nc)/(nc-1))
    sse.y <- sy^2/ny
    
    sse <- sse.x + sse.y
    sse2 <- sse^2
    
    ## test statistic
    stat <- (mx - my)/sqrt(sse)
    ## names(stat) <- rep("t", length(stat))
    
    ## approximate degrees of freedom (Welch-Satterthwaite)
    df <- sse2 / (sse.x^2 /(nx-1) + sse.y^2 / (ny-1))
    ## names(df) <- "df"
    
    ## p-value
    pval <- switch(alternative,
                   "two.sided" = 2*(1 - pt(abs(stat), df = df)),
                   "greater" = 1 - pt(stat, df = df),
                   "less" = pt(stat, df = df))
    
    list(statistic = stat,
                parameter = df,
                p.value = pval)
}




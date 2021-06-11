#' Welch T-tests for rows of a matrix
#'
#' @param X A \code{m x n} numeric matrix whose rows correspond to variables
#'   and columns to observations
#'
#' @param categ Either a numeric vector of \code{n} categories in \eqn{0, 1} for
#'   the observations, or a \code{n x B} matrix stacking \code{B} such vectors
#'   (typically permutations of an original vector of size \code{n})
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
#'   \item{estimate}{the mean difference between groups}}
#'   Each of these elements is a matrix of size \code{m x B}, coerced to a vector of length \code{m} if \code{B=1}
#'   
#' @author Pierre Neuvial
#' 
#' @details This function performs \code{m x B} Welch T tests on
#'   \code{n} observations using matrix operations. Its time complexity is
#'   \code{O(mBn)}. The code is much faster than using loops of 'apply'
#'   functions, especially for high-dimensional problems (small n and large m)
#'   because the overhead of the call to the 't.test' function is avoided and
#'   the code is vectorized
#'
#' @references B. L. Welch (1951), On the comparison of several mean values: an
#'   alternative approach. Biometrika, *38*, 330-336
#'
#' @export
#'
#' @examples
#'
#' m <- 300
#' n <- 38
#' mat <- matrix(rnorm(m*n), ncol=n)
#' categ <- rep(c(0, 1), times=c(27, n-27))
#' system.time(fwt <- rowWelchTests(mat, categ, alternative = "greater"))
#' str(fwt)
#'
#' # compare with ordinary t.test:
#' system.time(pwt <- apply(mat, 1, FUN=function(x) {
#'    t.test(x[categ==1], x[categ==0], alternative = "greater")$p.value
#' }))
#' all(abs(fwt$p.value-pwt) < 1e-10)  ## same results
#' 
#' # with several contrasts/permutations
#' B <- 100
#' categ_perm <- replicate(B, sample(categ))
#' system.time(fwt_perm <- rowWelchTests(mat, categ_perm, alternative = "greater"))
#' str(fwt_perm)
#' 
rowWelchTests <- function(X, categ, 
                          alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    stopifnot(all(categ %in% c(0, 1)))
    
    categ0 <- as.matrix(categ)
    categ1 <- 1 - categ0
    
    nX <- as.matrix(colSums(categ0))
    nY <- as.matrix(colSums(categ1))
    
    sumX <- X %*% categ0
    sumY <- X %*% categ1
    XX <- X * X
    sum2X <- XX %*% categ0
    sum2Y <- XX %*% categ1
    
    mX <- divide_cols(sumX, nX)
    mY <- divide_cols(sumY, nY)
    
    num <- sum2X - divide_cols(sumX^2, nX)
    sX <- sqrt(divide_cols(num, nX - 1))
    
    num <- sum2Y - divide_cols(sumY^2, nY)
    sY <- sqrt(divide_cols(num, nY - 1))
    
    suffWelchTests(mX, mY, sX, sY, nX, nY, alternative = alternative)
}

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
#' @return A list with elements:
#'   \describe{ 
#'   \item{statistic}{the value of the t-statistic}
#'   \item{parameter}{the degrees of freedom for the t-statistic}
#'   \item{p.value}{the p-value for the test}
#'   }
#' @author Pierre Neuvial
#'
#' @details In accordance with the implementation of \code{\link[stats]{t.test}} and friends, we
#'   follow the rule that 'alternative = "greater"' is the alternative that 'x'
#'   has a larger  mean than 'y'.
#'
#' @noRd
#' @importFrom stats pt
#' @examples
#'
#' # Reproducing the results of the 't.test' function
#'
#' x <- rnorm(100)
#' y <- rnorm(34, mean = 1)
#' target <- t.test(x, y)
#' swt <- suffWelchTests(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p.value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)
#' 
#' target <- t.test(x, y, alternative = "greater")
#' swt <- suffWelchTests(mean(x), mean(y), sd(x), sd(y), 
#'   length(x), length(y), alternative = "greater")
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p.value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)
#' 
suffWelchTests <- function(mx, my, sx, sy, nx, ny, 
                          alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    
    
    if (is.vector(mx)) {
        if (!all(is.vector(my), is.vector(sx), is.vector(sy),
                      is.vector(nx), is.vector(ny))) {
            stop("All numeric inputs should be of the same type (either vector or matrix)")
        }
        mx <- as.matrix(mx)
        my <- as.matrix(my)
        sx <- as.matrix(sx)
        sy <- as.matrix(sy)
        nx <- as.matrix(nx)
        ny <- as.matrix(ny)
    } else if (is.matrix(mx)) {
        if (!all(is.matrix(my), is.matrix(sx), is.matrix(sy),
                 is.matrix(nx), is.matrix(ny))) {
            stop("All numeric inputs should be of the same type (either vector or matrix)")
        }
    } else {
        stop("Numeric inputs should be vectors or matrices")
    }
    sse.x <- divide_cols(sx^2, nx)  ## sc <-> sqrt((sum2c - sumc^2/nc)/(nc-1))
    sse.y <- divide_cols(sy^2, ny)
    
    sse <- sse.x + sse.y
    sse2 <- sse^2
    
    ## test statistic
    stat <- (mx - my)/sqrt(sse)
    ## names(stat) <- rep("t", length(stat))
    
    ## approximate degrees of freedom (Welch-Satterthwaite)
    deno <- divide_cols(sse.x^2, nx-1) + divide_cols(sse.y^2, ny-1)
    df <- sse2/deno
    ## names(df) <- "df"
    
    # mean difference
    meanDiff <- mx - my
    
    # coerce to vector if single column
    if (ncol(stat) == 1) {
        stat <- stat[, 1]
        stopifnot(ncol(df) == 1)
        df <- df[, 1]
        stopifnot(ncol(meanDiff) == 1)
        meanDiff <- meanDiff[, 1]
    }
    
    ## p-value
    pval <- switch(alternative,
                   "two.sided" = 2*(1 - pt(abs(stat), df = df)),
                   "greater" = 1 - pt(stat, df = df),
                   "less" = pt(stat, df = df))
    
    list(statistic = stat,
         parameter = df,
         p.value = pval, 
         estimate = meanDiff)
}

divide_cols <- function(a, b) sweep(a, 2, b, "/")

categCheck <- function(categ, n) {
    name <- as.character(substitute(categ))
    if (length(categ) != n) {
        stop(name, " should be of length ", n, ", not", length(categ))
    }
    categ <- as.factor(categ)
    cats <- levels(categ)
    if (!identical(cats, c("0", "1"))) {
        stop("'", name, "' should consist only of '0' and '1'!")
    }
}

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
#' @importFrom stats pt
#' @noRd
#' @examples
#' # Reproducing the results of the 't.test' function
#'
#' x <- rnorm(100)
#' y <- rnorm(34, mean = 1)
#' target <- t.test(x, y)
#' swt <- suffWelchTests1(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p.value - target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter-target$parameter, check.attributes = FALSE)
#' 
#' target <- t.test(x, y, alternative = "greater")
#' swt <- suffWelchTests1(mean(x), mean(y), sd(x), sd(y), 
#'   length(x), length(y), alternative = "greater")
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p.value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)

#' @noRd
suffWelchTests1 <- function(mx, my, sx, sy, nx, ny, 
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
         p.value = pval,
         estimate = mx - my)
}

rowWelchTests1 <- function(mat, categ, alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    categCheck(categ, ncol(mat))
    
    sstats <- getSummaryStats(mat, categ = categ)
    Y <- sstats[["0"]]
    X <- sstats[["1"]]  ## as per the doc of t.test:
    ## 'alternative = "greater"' is the alternative that 'x' has a larger  mean
    ## than 'y'.
    swt <- suffWelchTests1(X[["mean"]], Y[["mean"]],
                          X[["sd"]], Y[["sd"]],
                          X[["n"]], Y[["n"]],
                          alternative = alternative)
    swt
}

#' Convert a matrix of observations into summary statistics
#' 
#' Convert a matrix of observations labelled into categories into summary 
#' statistics for each category
#' 
#' The following statistics are calculated: sums, sums of squares, means, 
#' standard deviations, sample sizes
#' 
#' @param mat A \code{m x n} numeric matrix whose rows correspond to variables
#'   and columns to observations
#'   
#' @param categ A vector of \code{n} categories for the observations
#'   
#' @return A list of \code{n} elements containing the above-described
#'   summary statistics for each category
#' 
#' @author Pierre Neuvial
#'   
#' @noRd
#' @examples
#' 
#' mat <- matrix(rnorm(3051*38), ncol = 38)
#' cls <- rep(c(0, 1), times = c(27, 11))
#' stats <- getSummaryStats(mat, categ = cls) 
getSummaryStats <- function(mat, categ) {
    stopifnot(ncol(mat) == length(categ))
    cats <- sort(unique(categ))
    res <- list()
    for (cc in seq(along = cats)) {
        ww <- which(categ == cats[cc])
        matc <- mat[, ww, drop = FALSE]
        
        sumc <- rowSums(matc)
        sum2c <- rowSums(matc^2)
        nc <- length(ww)
        mc <- sumc/nc
        sc <- sqrt((sum2c - sumc^2/nc)/(nc-1))
        
        res[[cc]] <- list(sum = sumc,
                          sum2 = sum2c,
                          n = nc,
                          mean = mc,
                          sd = sc)
    }
    names(res) <- cats
    return(res)
}

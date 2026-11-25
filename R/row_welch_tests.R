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
#'   \item{p_value}{the p-values for the tests}
#'   \item{estimate}{the mean difference between groups}}
#'   Each of these elements is a matrix of size \code{m x B}, coerced to a
#'   vector of length \code{m} if \code{B=1}
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
#' mat <- matrix(rnorm(m * n), ncol = n)
#' categ <- rep(c(0, 1), times = c(27, n - 27))
#' system.time(fwt <- row_welch_tests(mat, categ, alternative = "greater"))
#' str(fwt)
#'
#' # compare with ordinary t.test:
#' system.time(pwt <- apply(mat, 1, FUN = function(x) {
#'   t.test(x[categ == 1], x[categ == 0], alternative = "greater")$p.value
#' }))
#' all(abs(fwt$p_value - pwt) < 1e-10) ## same results
#'
#' # with several contrasts/permutations
#' B <- 100
#' categ_perm <- replicate(B, sample(categ))
#' system.time(fwt_perm <- row_welch_tests(mat, categ_perm,
#'   alternative = "greater"
#' ))
#' str(fwt_perm)
#'
row_welch_tests <- function(X, categ,
                            alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  stopifnot(all(categ %in% c(0, 1)))

  categ0 <- as.matrix(categ)
  categ1 <- 1 - categ0

  n_x <- as.matrix(colSums(categ0))
  n_y <- as.matrix(colSums(categ1))

  XX <- X * X
  sum_x <- X %*% categ0
  sum2_x <- (XX) %*% categ0
  m_x <- divide_cols(sum_x, n_x)

  num <- sum2_x - divide_cols(sum_x^2, n_x)
  rm(sum_x)
  rm(sum2_x)

  num <- abs(num) # fix rare corner cases where num ~= -1e-15
  s_x <- sqrt(divide_cols(num, n_x - 1))

  sum_y <- X %*% categ1
  sum2_y <- (XX) %*% categ1
  m_y <- divide_cols(sum_y, n_y)
  rm(XX)

  num <- sum2_y - divide_cols(sum_y^2, n_y)
  num <- abs(num) # fix rare corner cases where num ~= -1e-15
  rm(sum_y)
  rm(sum2_y)
  s_y <- sqrt(divide_cols(num, n_y - 1))
  rm(num)

  suff_welch_tests(m_x, m_y, s_x, s_y, n_x, n_y, alternative = alternative)
}

#' Welch test from sufficient statistics
#'
#' @param m_x A numeric value or vector, the sample average for condition "x"
#' @param m_y A numeric value or vector of the same length as 'm_x', the sample
#'   average for condition "y"
#' @param s_x A numeric value or vector of the same length as 'm_x', the
#'   standard deviation for condition "x"
#' @param s_y A numeric value or vector of the same length as 'm_x', the
#'   standard deviation for condition "y"
#' @param n_x A numeric value or vector of the same length as 'm_x', the sample
#'   size for condition "x"
#' @param n_y A numeric value or vector of the same length as 'm_x', the sample
#'   size for condition "y"
#' @param alternative A character string specifying the alternative hypothesis.
#' @return A list with elements:
#'   \describe{
#'   \item{statistic}{the value of the t-statistic}
#'   \item{parameter}{the degrees of freedom for the t-statistic}
#'   \item{p_value}{the p-value for the test}
#'   }
#' @author Pierre Neuvial
#'
#' @details In accordance with the implementation of \code{\link[stats]{t.test}}
#'   and friends, we follow the rule that 'alternative = "greater"' is the
#'   alternative that 'x' has a larger  mean than 'y'.
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
#' swt <- suff_welch_tests(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p_value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)
#'
#' target <- t.test(x, y, alternative = "greater")
#' swt <- suff_welch_tests(mean(x), mean(y), sd(x), sd(y),
#'   length(x), length(y),
#'   alternative = "greater"
#' )
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p_value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)
#'
suff_welch_tests <- function(m_x, m_y, s_x, s_y, n_x, n_y,
                             alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)


  if (is.vector(m_x)) {
    if (!all(
      is.vector(m_y), is.vector(s_x), is.vector(s_y),
      is.vector(n_x), is.vector(n_y)
    )) {
      stop(paste0("All numeric inputs should be of the same type",
        " (either vector or matrix)",
        sep = ""
      ))
    }
    m_x <- as.matrix(m_x)
    m_y <- as.matrix(m_y)
    s_x <- as.matrix(s_x)
    s_y <- as.matrix(s_y)
    n_x <- as.matrix(n_x)
    n_y <- as.matrix(n_y)
  } else if (is.matrix(m_x)) {
    if (!all(
      is.matrix(m_y), is.matrix(s_x), is.matrix(s_y),
      is.matrix(n_x), is.matrix(n_y)
    )) {
      stop(paste0("All numeric inputs should be of the same type",
        "(either vector or matrix)",
        sep = " "
      ))
    }
  } else {
    stop("Numeric inputs should be vectors or matrices")
  }
  sse_x <- divide_cols(s_x^2, n_x) ## sc <-> sqrt((sum2c - sumc^2/nc)/(nc-1))
  sse_y <- divide_cols(s_y^2, n_y)

  sse <- sse_x + sse_y
  sse2 <- sse^2

  ## test statistic
  statistics <- (m_x - m_y) / sqrt(sse)

  ## approximate degrees of freedom (Welch-Satterthwaite)
  deno <- divide_cols(sse_x^2, n_x - 1) + divide_cols(sse_y^2, n_y - 1)
  df <- sse2 / deno

  # mean difference
  mean_diff <- m_x - m_y

  # coerce to vector if single column
  if (ncol(statistics) == 1) {
    statistics <- statistics[, 1]
    stopifnot(ncol(df) == 1)
    df <- df[, 1]
    stopifnot(ncol(mean_diff) == 1)
    mean_diff <- mean_diff[, 1]
  }

  ## p-value
  p_values <- switch(alternative,
    "two.sided" = 2 * (1 - pt(abs(statistics), df = df)),
    "greater" = 1 - pt(statistics, df = df),
    "less" = pt(statistics, df = df)
  )

  list(
    statistic = statistics,
    parameter = df,
    p_value = p_values,
    estimate = mean_diff
  )
}

divide_cols <- function(a, b) sweep(a, 2, b, "/")

#' Welch test from sufficient statistics
#'
#' @param m_x A numeric value or vector, the sample average for condition "x"
#' @param m_y A numeric value or vector of the same length as 'm_x', the sample
#'   average for condition "y"
#' @param s_x A numeric value or vector of the same length as 'm_x', the
#'   standard deviation for condition "x"
#' @param s_y A numeric value or vector of the same length as 'm_x', the
#'   standard deviation for condition "y"
#' @param n_x A numeric value or vector of the same length as 'm_x', the sample
#'   size for condition "x"
#' @param n_y A numeric value or vector of the same length as 'm_x', the sample
#'   size for condition "y"
#' @param alternative A character string specifying the alternative hypothesis.
#'   Currently only "two.sided" is implemented
#' @return A list with elements:
#'   \describe{ 
#'   \item{statistic}{the value of the t-statistic}
#'   \item{p_value}{the p-value for the test}
#'   \item{parameter}{the degrees of freedom for the t-statistic}
#'   }
#' @author Pierre Neuvial
#'
#' @details In accordance with the implementation of \code{\link[stats]{t.test}}
#'   and friends, we follow the rule that 'alternative = "greater"' is the
#'   alternative that 'x' has a larger  mean than 'y'.
#'
#' @importFrom stats pt
#' @noRd
#' @examples
#' # Reproducing the results of the 't.test' function
#'
#' x <- rnorm(100)
#' y <- rnorm(34, mean = 1)
#' target <- t.test(x, y)
#' swt <- suff_welch_tests1(mean(x), mean(y),
#'                          sd(x), sd(y), length(x), length(y))
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p_value - target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter - target$parameter, check.attributes = FALSE)
#'
#' target <- t.test(x, y, alternative = "greater")
#' swt <- suff_welch_tests1(mean(x), mean(y), sd(x), sd(y),
#'   length(x), length(y),
#'   alternative = "greater"
#' )
#' all.equal(swt$statistic, target$statistic, check.attributes = FALSE)
#' all.equal(swt$p_value, target$p.value, check.attributes = FALSE)
#' all.equal(swt$parameter, target$parameter, check.attributes = FALSE)
#'
#' @noRd
suff_welch_tests1 <- function(m_x, m_y, s_x, s_y, n_x, n_y,
                              alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  #    if (alternative != "two.sided") stop("altenative", alternative,
  #                                         "not implemented yet!")

  ## sanity checks
  p <- length(m_x)
  stopifnot(length(m_y) == p)
  stopifnot(length(s_x) == p)
  stopifnot(length(s_y) == p)
  stopifnot(length(n_x) %in% c(1, p))
  stopifnot(length(n_y) %in% c(1, p))

  sse_x <- s_x^2 / n_x
  sse_y <- s_y^2 / n_y

  sse <- sse_x + sse_y
  sse2 <- sse^2

  ## test statistic
  statistics <- (m_x - m_y) / sqrt(sse)

  ## approximate degrees of freedom (Welch-Satterthwaite)
  df <- sse2 / (sse_x^2 / (n_x - 1) + sse_y^2 / (n_y - 1))

  ## p-value
  p_value <- switch(alternative,
    "two.sided" = 2 * (1 - pt(abs(statistics), df = df)),
    "greater" = 1 - pt(statistics, df = df),
    "less" = pt(statistics, df = df)
  )

  list(
    statistic = statistics,
    parameter = df,
    p_value = p_value,
    estimate = m_x - m_y
  )
}

row_welch_tests1 <- function(mat, categ,
                             alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  categ_check(categ, ncol(mat))

  summary_statistics <- summary_statistics(mat, categ = categ)
  Y <- summary_statistics[["0"]]
  X <- summary_statistics[["1"]] ## as per the doc of t.test:
  ## 'alternative = "greater"' is the alternative that 'x' has a larger  mean
  ## than 'y'.
  swt <- suff_welch_tests1(X[["mean"]], Y[["mean"]],
    X[["sd"]], Y[["sd"]],
    X[["n"]], Y[["n"]],
    alternative = alternative
  )
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
#' mat <- matrix(rnorm(3051 * 38), ncol = 38)
#' cls <- rep(c(0, 1), times = c(27, 11))
#' stats <- summary_statistics(mat, categ = cls)
summary_statistics <- function(mat, categ) {
  stopifnot(ncol(mat) == length(categ))
  cats <- sort(unique(categ))
  res <- list()
  for (cc in seq(along = cats)) {
    ww <- which(categ == cats[cc])
    matc <- mat[, ww, drop = FALSE]

    sumc <- rowSums(matc)
    sum2c <- rowSums(matc^2)
    nc <- length(ww)
    mc <- sumc / nc
    sc <- sqrt((sum2c - sumc^2 / nc) / (nc - 1))

    res[[cc]] <- list(
      sum = sumc,
      sum2 = sum2c,
      n = nc,
      mean = mc,
      sd = sc
    )
  }
  names(res) <- cats
  return(res)
}

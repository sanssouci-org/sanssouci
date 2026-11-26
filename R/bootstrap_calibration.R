#' Mass-univariate inference for contrasts in a linear model
#'
#' Compute the marginal t-statistics for a set of contrasts and their
#' (two-sided) p-value.

#' @param Y A data matrix of size $n$ observations (in row) and $D$ features in
#'   columns
#' @param X A design matrix of size $n$ observations (in row) and $p$ variables
#'   (in columns)
#' @param C A contrast matrix of size $L$ tested contrasts (in row) and $p$
#'   columns corresponding to the parameters to be tested
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#'
#' @return A list with elements:
#' \describe{
#'   \item{statistic}{A \eqn{L \times D} matrix of test statistics}
#'   \item{p_value}{A \eqn{L \times D} matrix of p-values}
#'   \item{estimate}{A \eqn{n \times D} matrix of parameter estimates}
#'   \item{residual}{A \eqn{n \times D} matrix of residuals}
#' }
#' @details Based on a python implementation available in the \code{pyperm}
#'   package: \url{https://github.com/sjdavenport/pyperm/}
#'
#' @references Davenport, S., Thirion, B., & Neuvial, P. (2025). FDP control in
#'   mass-univariate linear models using the residual bootstrap. Electronic
#'   Journal of Statistics, 19(1), 1313-1336.
#'
#' @export
#'
#' @examples
#' N <- 100
#' p <- 2
#' D <- 2
#' X <- matrix(0, nrow = N, ncol = p)
#' X[, 1] <- 1
#' X[, -1] <- runif(N * (p - 1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N * D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' res_lm <- lm_test(Y = Y, X = X, C = C)
lm_test <- function(Y, X, C, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  .check_lm_test(Y, X, C)

  df <- nrow(Y) - qr(X)$rank

  ## estimation of model parameters
  xtx <- crossprod(X)
  xty <- crossprod(X, Y)
  xtx_inv <- solve(xtx)
  estimates <- crossprod(xtx_inv, xty)

  ## residuals
  residuals <- Y - X %*% estimates

  ## estimation of residual standard deviation
  sigma_est <- sqrt(colSums(residuals^2) / df)

  ## test statistic
  cxtx_invct <- C %*% xtx_inv %*% t(C) # matrix L x L
  cxtx_diag <- sqrt(diag(cxtx_invct)) # vector L
  denom <- outer(cxtx_diag, sigma_est, "*") # matrix L x D
  statistics <- C %*% estimates / denom # matrix (L x D)

  ## p-values
  p_values <- switch(alternative,
    "two.sided" = 2 * (1 - pt(abs(statistics), df = df)),
    "greater" = 1 - pt(statistics, df = df),
    "less" = pt(statistics, df = df)
  )

  return(list(
    p_value = p_values,
    statistic = statistics,
    estimate = estimates,
    residual = residuals
  ))
}

#' Mass-univariate bootstrap-based inference for contrasts in a linear model
#'
#' Compute the marginal null t-statistics for a set of contrasts and their
#' (two-sided) p-value by bootstrapping the residuals
#'
#' @inheritParams lm_test
#' @param B An integer value, the number of bootstraps
#' @param replace A Boolean value. If TRUE (default) then the residuals are
#'   sampled with replacement (i.e. a bootstrap), if FALSE then they are sampled
#'   without replacement resulting in a permutation of the data
#'
#' @return An array of permuted p-values of dimensions \eqn{D \times L \times B}
#' @details Performs \code{lm_test} for each permutation. Based on a python
#'   implementation available in the \code{pyperm} package:
#'   \url{https://github.com/sjdavenport/pyperm/}
#'
#' @references Davenport, S., Thirion, B., & Neuvial, P. (2025). FDP control in
#'   mass-univariate linear models using the residual bootstrap. Electronic
#'   Journal of Statistics, 19(1), 1313-1336.
#'
#' @export
#'
#' @examples
#' N <- 100
#' p <- 2
#' D <- 2
#' X <- matrix(0, nrow = N, ncol = p)
#' X[, 1] <- 1
#' X[, -1] <- runif(N * (p - 1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N * D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' res_lm <- bootstrap_permutation(Y = Y, X = X, C = C, B = 10)
bootstrap_permutation <- function(Y, X, C,
                                  alternative = c(
                                    "two.sided", "less",
                                    "greater"
                                  ), B = 1000,
                                  replace = TRUE) {
  .check_lm_test(Y, X, C)
  if (!is.numeric(B) || B < 0) {
    stop("B must be a positive numeric value",
      call. = FALSE
    )
  }
  alternative <- match.arg(alternative)

  n <- nrow(Y)
  D <- ncol(Y)
  L <- nrow(C)

  ## perform LM test to obtain esiduals
  res_lm <- lm_test(Y = Y, X = X, C = C, alternative = alternative)
  epsilon_hat <- res_lm$residual

  ## Bootstrapping of residuals
  epsilon_perm <- array(NA, dim = c(dim(epsilon_hat), B))
  for (b in 1:B) {
    shuffle_idx <- sample.int(n, replace = replace)
    epsilon_perm[, , b] <- epsilon_hat[shuffle_idx, ]
  }

  ## computation of permuted expression Y_n^{(b)}
  y_perm <- epsilon_perm

  ## Test contrasts for each permuted expression matrix
  perm_p_values <- array(NA, dim = c(D, L, B))
  for (b in 1:B) {
    perm_p_values[, , b] <- lm_test(
      Y = as.matrix(y_perm[, , b]),
      X = X, C = C,
      alternative = alternative
    )$p_value
  }

  return(perm_p_values)
}

#' Mass-univariate bootstrap-based inference for contrasts in a linear model
#'
#' Compute the marginal null t-statistics for a set of contrasts and their
#' (two-sided) p-value by bootstrapping the residuals
#'
#' @inheritParams lm_test
#' @param alternative  A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#' @param groups A numeric matrix of \eqn{n} rows and \eqn{B} columns values in
#' \eqn{1, ..., n},  the indicator
#' of the sample used in the test.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{p_value}{A vector of p-values}
#'   \item{statistic}{A vector of test statistic}
#'   \item{estimate}{A vector of OLS estimators}
#'   }
#'   Each of these elements is a matrix of size \code{m x B},
#'   coerced to a vector of length \code{m} if \code{B=1}
#' @export
row_lm_test <- function(Y, X, C,
                        alternative = c("two.sided", "less", "greater"),
                        groups = matrix(seq_len(ncol(Y)), ncol = 1)) {
  Y <- t(Y)
  .check_lm_test(Y, X, C)
  if (nrow(X) != nrow(groups)) {
    stop("groups must be a matrix with $n = nrow(X)$ rows",
      call. = FALSE
    )
  }
  alternative <- match.arg(alternative)

  n <- nrow(Y)
  D <- ncol(Y)
  L <- nrow(C)
  B <- ncol(groups)

  res_lm <- lm_test(Y = Y, X = X, C = C, alternative = alternative)
  residuals <- res_lm$residual

  if (all(groups == 1:n) && ncol(matrix(groups)) == 1) {
    estimate <- C %*% res_lm$estimate
    return(list(
      p_value = res_lm$p_value,
      statistic = res_lm$statistic,
      estimate = estimate
    ))
  }

  ## Bootstrapping of residuals
  perm_statistics <- array(NA, dim = c(D, L, B))
  perm_p_values <- array(NA, dim = c(D, L, B))
  perm_estimates <- array(NA, dim = c(D, L, B))
  for (b in 1:B) {
    shuffle_idx <- groups[, b]
    y_perm <- residuals[shuffle_idx, ]
    res_perm <- lm_test(
      Y = as.matrix(y_perm), X = X, C = C,
      alternative = alternative
    )
    perm_statistics[, , b] <- res_perm$statistic
    perm_p_values[, , b] <- res_perm$p_value
    perm_estimates[, , b] <- C %*% res_perm$estimate
  }
  ## transform 3 dimensional problem (D, L, B) into a matrix (D*L, B)
  perm_p_values_ <- matrix(perm_p_values, nrow = D * L, ncol = B)
  perm_statistics_ <- matrix(perm_statistics, nrow = D * L, ncol = B)
  perm_estimates_ <- matrix(perm_estimates, nrow = D * L, ncol = B)

  return(list(
    statistic = perm_statistics_,
    p_value = perm_p_values_,
    estimate = perm_estimates_
  ))
}

#' Calibration of post hoc bound using bootstrap permutations
#'
#' Compute by bootstraping a Joint Error Rate controlling threshold family
#' associated to a set of contrast in a linear model.
#'
#' @inheritParams bootstrap_permutation
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#' @param family A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes", or "Oracle". "Linear" and "Simes"
#'   families are identical.
#'
#'   - Simes/Linear: The classical family of thresholds introduced by
#'   Simes (1986). This family yields JER control if the test statistics are
#'   positively dependent (PRDS) under H0.
#'
#'   - Beta: A family of thresholds that achieves marginal kFWER control under
#'   independence
#'
#'   - Oracle A family such that the associated bounds correspond to the true
#'   numbers/proportions of true/false positives. "truth" must be available in
#'   object$input$truth.
#'
#' @return A list with elements:
#' \describe{
#'   \item{thr}{A numeric vector of length K, such that the estimated
#'   probability that there exists an index k between 1 and K such that the k-th
#'   maximum of the test statistics of is greater than `thr[k]`, is less than
#'   alpha}
#'   \item{piv_stat}{A vector of `B` pivotal statitics}
#'   \item{lambda}{A numeric value, the result of the calibration}
#' }
#'
#' @references Davenport, S., Thirion, B., & Neuvial, P. (2025). FDP control in
#'   mass-univariate linear models using the residual bootstrap. Electronic
#'   Journal of Statistics, 19(1), 1313-1336.
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc
#'   confidence bounds on false positives using reference families.
#'
#' @export
#'
#' @examples
#' N <- 100
#' p <- 2
#' D <- 2
#' X <- matrix(0, nrow = N, ncol = p)
#' X[, 1] <- 1
#' X[, -1] <- runif(N * (p - 1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N * D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' res_lm <- calibration_bootstap(Y = Y, X = X, C = C, B = 10)
calibration_bootstap <- function(Y, X, C,
                                 alternative = c(
                                   "two.sided", "less", "greater"
                                 ),
                                 B = 1000, alpha = 0.05,
                                 family = c(
                                   "Simes", "Linear", "Beta", "Oracle"
                                 )) {
  .check_lm_test(Y, X, C)
  if (!is.numeric(B) || B < 0) {
    stop("B must be a positive numeric value",
      call. = FALSE
    )
  }
  family <- match.arg(family)
  alternative <- match.arg(alternative)

  D <- dim(Y)[2]
  L <- dim(C)[1]

  ## perform permuted pvalues
  perm_p_values <- bootstrap_permutation(
    Y = Y, X = X, C = C, B = B, replace = TRUE,
    alternative = alternative
  )

  ## transform 3 dimensional problem (D, L, B) into a matrix (D*L, B)
  perm_p_values_ <- matrix(perm_p_values, nrow = D * L, ncol = B)

  ## The next steps of the calibration are already implemented in sanssouci
  cal <- calibrate0(perm_p_values_, m = D * L, alpha = alpha, family = family)
  return(cal)
}

.check_lm_test <- function(Y, X, C) {
  vars <- list(Y, X, C)
  names <- as.character(substitute(list(Y, X, C)))[-1]
  for (i in seq_along(vars)) {
    if (!is.matrix(vars[[i]])) {
      stop(sprintf("'%s' must be a matrix", names[i]), call. = FALSE)
    }
  }

  if (nrow(Y) != nrow(X)) {
    stop("nrow(Y) must be equal to nrow(X). This is
                              corresponding to $n$ the number of observations",
      call. = FALSE
    )
  }
  if (ncol(X) != ncol(C)) {
    stop("ncol(X) must be equal to ncol(C). This is
                              corresponding to $p$ the number of variables",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

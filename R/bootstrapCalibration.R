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
#'   \item{epsilon_est}{A \eqn{n \times D} matrix of residuals}
#'   \item{stat_test}{A \eqn{L \times D} matrix of test statistics}
#'   \item{p.value}{A \eqn{L \times D} matrix of p-values}
#'   \item{beta_est}{A \eqn{n \times D} matrix of parameter estimates}
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
#' X[, -1] <- runif(N*(p-1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' resLM <- lm_test(Y = Y, X = X, C = C)
lm_test <- function(Y, X, C, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  .check_lm_test(Y, X, C)
  
  df <- nrow(Y) - qr(X)$rank
  
  ## estimation of model parameters
  XtX <- crossprod(X) # t(X) %*% X
  XtY <- crossprod(X, Y) # t(X) %*% Y
  XtX_inv <- solve(XtX)
  beta_est <- crossprod(XtX_inv, XtY)
  
  ## residuals
  epsilon_est <- Y - X %*% beta_est
  
  ## estimation of residual standard deviation
  sigma_est <- sqrt(colSums(epsilon_est^2) / df)
  
  ## test statistic
  CXtX_invCt <- C %*% XtX_inv %*% t(C)      # matrix L x L
  CXtX_diag <- sqrt(diag(CXtX_invCt))       # vector L
  denom <- outer(CXtX_diag, sigma_est, "*") # matrix L x D
  stat_test <- C %*% beta_est / denom       # matrix (L x D)
  
  ## p-values
  pval <- switch(alternative,
                 "two.sided" = 2 * (1 - pt(abs(stat_test), df = df)),
                 "greater" = 1 - pt(stat_test, df = df),
                 "less" = pt(stat_test, df = df)
  )
  
  return(list(
    epsilon_est = epsilon_est,
    stat_test = stat_test,
    p.value = pval,
    beta_est = beta_est
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
#' X[, -1] <- runif(N*(p-1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N * D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' resLM <- bootstrap_permutation(Y = Y, X = X, C = C, B = 10)
bootstrap_permutation <- function(Y, X, C, 
                                  alternative = c("two.sided", "less",
                                                  "greater"), B = 1000, 
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
  resLM <- lm_test(Y = Y, X = X, C = C, alternative = alternative)
  epsilon_hat <- resLM$epsilon_est
  
  ## Bootstrapping of residuals
  epsilon_perm <- array(NA, dim = c(dim(epsilon_hat), B))
  for (b in 1:B) {
    shuffle_idx <- sample.int(n, replace = replace)
    epsilon_perm[, , b] <- epsilon_hat[shuffle_idx, ]
  }
  # Xhatbeta <- X %*% resLM$beta_est
  
  ## computation of permuted expression Y_n^{(b)}
  Y_perm <- epsilon_perm
  # Y_perm <- sweep(epsilon_perm, 1:2, Xhatbeta, "+")
  
  ## Test contrasts for each permuted expression matrix
  pval_perm <- array(NA, dim = c(D, L, B))
  for (b in 1:B) {
    pval_perm[, , b] <- lm_test(
      Y = as.matrix(Y_perm[, , b]),
      X = X, C = C, 
      alternative = alternative
    )$p.value
  }
  
  return(pval_perm)
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
#'   \item{statistic}{the value of the t-statistics}
#'   \item{p.value}{the p-values for the tests}
#'   Each of these elements is a matrix of size \code{m x B}, 
#'   coerced to a vector of length \code{m} if \code{B=1}
#' }
#' @export
row_lm_test <- function(Y, X, C,
                        alternative = c("two.sided", "less", "greater"), 
                        groups = matrix(1:ncol(Y), ncol = 1)){
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
  
  resLM <- lm_test(Y = Y, X = X, C = C, alternative = alternative)
  epsilon_hat <- resLM$epsilon_est
  
  if(all(groups == 1:n) & ncol(matrix(groups)) == 1){
    estimate <- C %*% resLM$beta_est
    return(list(p.value = resLM$p.value, statistic = resLM$stat_test, 
                estimate = estimate))
  }
  
  ## Bootstrapping of residuals
  pval_perm <- array(NA, dim = c(D, L, B))
  stat_perm <- array(NA, dim = c(D, L, B))
  estimate_perm <- array(NA, dim = c(D, L, B))
  for (b in 1:B) {
    shuffle_idx <- groups[,b]
    Y_perm <- epsilon_hat[shuffle_idx, ]
    res_perm <- lm_test(Y = as.matrix(Y_perm), X = X, C = C, 
                        alternative = alternative)
    pval_perm[, , b] <- res_perm$p.value
    stat_perm[, , b] <- res_perm$stat_test
    estimate_perm[, , b] <- C %*% res_perm$beta_est
  }
  ## transform 3 dimensional problem (D, L, B) into a matrix (D*L, B)
  pval_perm_matrix <- matrix(pval_perm, nrow = D * L, ncol = B)
  stat_perm_matrix <- matrix(stat_perm, nrow = D * L, ncol = B)
  estimate_perm_matrix <- matrix(estimate_perm, nrow = D * L, ncol = B)
  
  return(list(p.value = pval_perm_matrix, statistic = stat_perm_matrix, 
              estimate = estimate_perm_matrix))
  
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
#'   - Simes/Linear: The classical family of thresholds introduced by Simes (1986).
#'   This family yields JER control if the test statistics are positively
#'   dependent (PRDS) under H0.
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
#'   \item{thr}{A numeric vector of length K, such that the estimated probability that
#' there exists an index k between 1 and K such that the k-th maximum of the
#' test statistics of is greater than `thr[k]`, is less than alpha}
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
#' X <- matrix(0,nrow = N, ncol = p)
#' X[, 1] <- 1
#' X[, -1] <- runif(N*(p-1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' resLM <- calibration_bootstap(Y = Y, X = X, C = C, B = 10)
calibration_bootstap <- function(Y, X, C, 
                                 alternative = c("two.sided", "less", "greater"), 
                                 B = 1000, alpha = 0.05, 
                                 family = c("Simes", "Linear", "Beta", "Oracle")) {
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
  pval_perm <- bootstrap_permutation(Y = Y, X = X, C = C, B = B, replace = TRUE, 
                                     alternative = alternative) 
  
  ## transform 3 dimensional problem (D, L, B) into a matrix (D*L, B)
  pval_perm_matrix <- matrix(pval_perm, nrow = D * L, ncol = B)
  
  ## The next steps of the calibration are already implemented in sanssouci
  return(calibrate0(pval_perm_matrix,
                    m = D * L, alpha = alpha,
                    family = family
  ))
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

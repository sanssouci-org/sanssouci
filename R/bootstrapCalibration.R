#' Perform statistical test from contrast matrix for Linear model.
#'
#' @param Y Data matrix of size $n$ observations (in row) and $D$
#' features in columns
#' @param X Design matrix of size $n$ observations (in row) and $p$ variables
#' (in columns)
#' @param C Contrast matrix of size $L$ tested contrasts (in row) and $p$ columns
#' corresponding to the parameters to be tested
#'
#' @return A list with elements:
#'     \describe{
#'   \item{epsilon_est}{\eqn{n \times D} matrix of residuals}
#'   \item{stat_test}{\eqn{L \times D} matrix of test statistics}
#'   \item{pvalues}{\eqn{L \times D} matrix of p-values}
#'   \item{beta_est}{\eqn{n \times D} matrix of parameter estimates}  
#' }
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
lm_test <- function(Y, X, C) {
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
  pvaleurs <- 2 * pt(-abs(stat_test), df = df)

  return(list(
    epsilon_est = epsilon_est,
    stat_test = stat_test,
    pvalues = pvaleurs,
    beta_est = beta_est
  ))
}

#' Permuted p-values using bootstrap for linear model
#' @param Y Expression matrix of size $n$ observations (in row) and $D$
#' dimensions/mesures in columns
#' @param X Design matrix of size $n$ observations (in row) and $p$ variables
#' (in columns)
#' @param C Contrast matrix of size $L$ tested contrast (in row) and $p$ columns
#' corresponding to te tested parameters
#' @param B scalar, number of bootstrap
#' @param replace Bool. If TRUE (default) then the residuals are sampled with
#' replacement (i.e. a bootstrap), if FALSE then they are sampled without
#' replacement resulting in a permutation of the data
#'
#' @return matrix of permuted $p$-values. It is an array of dimensions 
#' \eqn{D \times L \times B}.
#' @export
#'
#' @examples
#' N = 100
#' p = 2
#' D = 2
#' X <- matrix(0,nrow = N, ncol = p)
#' X[,1] <- 1
#' X[,-1] <- runif(N*(p-1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' resLM <- bootstrap_permutation(Y = Y, X = X, C = C, B = 10)
bootstrap_permutation <- function(Y, X, C, B = 1000, replace = TRUE) {
  set.seed(13012024)
  .check_lm_test(Y, X, C)
  if (!is.numeric(B) || B < 0) {
    stop("B must be a positive numeric value",
      call. = FALSE
    )
  }

  n <- nrow(Y)
  D <- ncol(Y)
  L <- nrow(C)

  ## perform LM test to obtain esiduals
  resLM <- lm_test(Y = Y, X = X, C = C)
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
      X = X, C = C
    )$pvalues
  }

  return(pval_perm)
}

#' Calibration of post hoc bound using bootstrap permutations
#'
#' @param Y Expression matrix of size $n$ observations (in row) and $D$
#' dimensions/mesures in columns
#' @param X Design matrix of size $n$ observations (in row) and $p$ variables
#' (in columns)
#' @param C Contrast matrix of size $L$ tested contrast (in row) and $p$ columns
#' corresponding to te tested parameters
#' @param B scalar, number of bootstrap
#'
#' @return A list with elements
#' - thr: A numeric vector of length K, such that the estimated probability that
#' there exists an index k between 1 and K such that the k-th maximum of the
#' test statistics of is greater than `thr[k]`, is less than alpha
#' - piv_stat: A vector of `B` pivotal statitics
#' - lambda: A numeric value, the result of the calibration
#' @export
#'
#' @examples
#' N = 100
#' p = 2
#' D = 2
#' X <- matrix(0,nrow = N, ncol = p)
#' X[,1] <- 1
#' X[,-1] <- runif(N*(p-1), min = 0, max = 3)
#' beta <- matrix(0, nrow = p, ncol = D)
#' epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
#' Y <- X %*% beta + epsilons
#' C <- diag(p)
#' resLM <- calibration_bootstap(Y = Y, X = X, C = C, B = 10)
calibration_bootstap <- function(Y, X, C, B = 1000) {
  .check_lm_test(Y, X, C)
  if (!is.numeric(B) || B < 0) {
    stop("B must be a positive numeric value",
      call. = FALSE
    )
  }

  D <- dim(Y)[2]
  L <- dim(C)[1]

  ## perform permuted pvalues
  pval_perm <- bootstrap_permutation(Y = Y, X = X, C = C, B = B, replace = TRUE)

  ## transform 3 dimensional problem (D, L, B) into a matrix (D*L, B)
  pval_perm_martix <- matrix(pval_perm, nrow = D * L, ncol = B)

  ## The next steps of the calibration are already implemented in sanssouci
  return(sanssouci::calibrate0(pval_perm_martix,
    m = D * L, alpha = 0.05,
    family = "Linear"
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

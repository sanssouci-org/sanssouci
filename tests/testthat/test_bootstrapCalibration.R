library(testthat)
library(reticulate)

test_that('Consistency between `sanssouci::lm_test` and stat::lm', {
  
  ## If X is continuous
  DD <- c(1, 10)
  PP <- c(2, 10)
  N <- 100
  for (D in DD) {
    for (p in PP) {
      X <- matrix(0, nrow = N, ncol = p)
      X[,1] <- 1
      X[, -1] <- runif(N*(p - 1), min = 0, max = 3)
      beta <- matrix(0, nrow = p, ncol = D)
      epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
      Y <- X %*% beta + epsilons
      C <- diag(p)
      
      resLM <- lm_test(Y = Y, X = X, C = C)
      
      expect_type(resLM, "list")
      expect_true("epsilon_est" %in% names(resLM))
      expect_true("stat_test" %in% names(resLM))
      expect_true("pvalues" %in% names(resLM))
      expect_true("beta_est" %in% names(resLM))
      expect_equal(dim(resLM$epsilon_est), c(nrow(Y), ncol(Y)))
      expect_equal(dim(resLM$stat_test), c(nrow(C), ncol(Y)))
      expect_equal(dim(resLM$pvalues), c(nrow(C), ncol(Y)))
      expect_equal(dim(resLM$beta_est), c(ncol(X), ncol(Y)))
      
      for (j in 1:ncol(Y)) {
        fit <- lm(Y[, j] ~ X[, -1])  
        s <- summary(fit)
        p_lm <- unname(coef(s)[, 4])  # p-value
        expect_equal(resLM$pvalues[, j], p_lm, tolerance = 1e-6)
        stat_lm <- unname(coef(s)[, 3])  # stat
        expect_equal(resLM$stat_test[, j], stat_lm, tolerance = 1e-6)
        beta_lm <- unname(coef(s)[, "Estimate"]) # estimated beta
        expect_equal(resLM$beta_est[, j], beta_lm, tolerance = 1e-6)
      }
    }
  }
  
  ## If X is qualitative
  N <- 100
  for (D in DD) {
    for (p in PP) {
      X <- matrix(0,nrow = N, ncol = p)
      X[, 1] <- 1
      X[, -1] <- sample(0:1, N*(p - 1), replace = TRUE)
      beta <- matrix(0, nrow = p, ncol = D)
      epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
      Y <- X %*% beta + epsilons
      C <- diag(p)
      
      resLM <- lm_test(Y = Y, X = X, C = C)
      
      expect_type(resLM, "list")
      expect_true("epsilon_est" %in% names(resLM))
      expect_true("stat_test" %in% names(resLM))
      expect_true("pvalues" %in% names(resLM))
      expect_true("beta_est" %in% names(resLM))
      expect_equal(dim(resLM$epsilon_est), c(nrow(Y), ncol(Y)))
      expect_equal(dim(resLM$stat_test), c(nrow(C), ncol(Y)))
      expect_equal(dim(resLM$pvalues), c(nrow(C), ncol(Y)))
      expect_equal(dim(resLM$beta_est), c(ncol(X), ncol(Y)))
      
      for (j in 1:ncol(Y)) {
        fit <- lm(Y[, j] ~ X[, -1])  
        s <- summary(fit)
        p_lm <- unname(coef(s)[, 4])  # p-value 
        expect_equal(resLM$pvalues[, j], p_lm, tolerance = 1e-6)
        stat_lm <- unname(coef(s)[, 3])  # stat
        expect_equal(resLM$stat_test[, j], stat_lm, tolerance = 1e-6)
        beta_lm <- unname(coef(s)[,"Estimate"]) # estimated beta
        expect_equal(resLM$beta_est[, j], beta_lm, tolerance = 1e-6)
      }
    }
  }
  
  
  ## Error class of object
  N <- 10
  P <- 2
  D <- 2
  L <- 1
  X <- matrix(0,nrow = N, ncol = p)
  Y <- matrix(1, nrow = N, ncol = D)
  C <- matrix(1, nrow = L, ncol = D)
  error_list <- list("a", c("a","b"), 1, 1:3, list("a","b"), array(0, dim = c(2,3,2)))
  for (error in error_list) {
    expect_error(lm_test(Y = error, X = X, C = C), 
                 regexp = "'Y' must be a matrix")
    expect_error(lm_test(Y = Y, X = error, C = C),
                 regexp = "'X' must be a matrix")
    expect_error(lm_test(Y = Y, X = X, C = error),
                 regexp = "'C' must be a matrix")
  }
  
  # test mismatch in number of observation
  X <- matrix(0,nrow = N, ncol = p)
  Y <- matrix(1, nrow = N+1, ncol = D)
  C <- matrix(1, nrow = L, ncol = D)
  expect_error(lm_test(Y = Y, X = X, C = C))
  
  # test mismatch in number of variables
  X <- matrix(0,nrow = N, ncol = p)
  Y <- matrix(1, nrow = N, ncol = D)
  C <- matrix(1, nrow = L, ncol = D+1)
  expect_error(lm_test(Y = Y, X = X, C = C))
})


test_that("`bootstrap_permutation`", {
  DD <- c(1, 10)
  PP <- c(2, 10)
  BB <- c(10, 100)
  N <- 100
  for (D in DD) {
    for (p in PP) {
      for (B in BB) {
        
        X <- matrix(0, nrow = N, ncol = p)
        X[,1] <- 1
        X[,-1] <- runif(N*(p - 1), min = 0, max = 3)
        beta <- matrix(0, nrow = p, ncol = D)
        epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
        Y <- X %*% beta + epsilons
        C <- diag(p)
        
        pval_perm <- bootstrap_permutation(Y = Y, X = X, C = C, B, 
                                           replace = TRUE )
        
        expect_equal(class(pval_perm), "array")
        expect_equal(dim(pval_perm), c(ncol(Y), nrow(C), B ))
        expect_true(all(pval_perm <= 1))
        expect_true(all(pval_perm >= 0))
      }
    }
  }
  
  ## Error class of object
  N <- 10
  P <- 2
  D <- 2
  L <- 1
  B <- 10
  X <- matrix(0, nrow = N, ncol = p)
  Y <- matrix(1, nrow = N, ncol = D)
  C <- matrix(1, nrow = L, ncol = D)
  error_list <- list("a", 
                     c("a","b"), 
                     1, 
                     1:3, 
                     list("a","b"), 
                     array(0, dim = c(2,3,2)))
  for (error in error_list) {
    expect_error(bootstrap_permutation(Y = error, X = X, C = C, B = B), 
                 regexp = "'Y' must be a matrix")
    expect_error(bootstrap_permutation(Y = Y, X = error, C = C, B = B),
                 regexp = "'X' must be a matrix")
    expect_error(bootstrap_permutation(Y = Y, X = X, C = error, B = B),
                 regexp = "'C' must be a matrix")
  }
  
  # test mismatch in number of observation
  X <- matrix(0, nrow = N, ncol = p)
  Y <- matrix(1, nrow = N + 1, ncol = D)
  C <- matrix(1, nrow = L, ncol = D)
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = B))
  
  # test mismatch in number of variables
  X <- matrix(0,nrow = N, ncol = p)
  Y <- matrix(1, nrow = N, ncol = D)
  C <- matrix(1, nrow = L, ncol = D + 1)
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = B))
  
  # test values of B
  X <- matrix(0,nrow = N, ncol = p)
  Y <- matrix(1, nrow = N, ncol = D)
  C <- matrix(1, nrow = L, ncol = D)
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = "a"))
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = matrix(1, 
                                                                     nrow = 2, 
                                                                     ncol = 3)))
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = list("a")))
  expect_error(bootstrap_permutation(Y = Y, X = X, C = C, B = -5))
})



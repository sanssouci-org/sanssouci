skip_on_ci()
skip_on_cran()
skip(message = "Skip test of R - Python consistency for bootstrap calibration")

library(reticulate)
if (!virtualenv_exists("env_test")) {
  virtualenv_create("env_test")
}
use_virtualenv("env_test", required = TRUE)

# Install pyperm
if (!py_module_available("pyperm")) {
  virtualenv_install("env_test", 
                     packages = c("numpy==1.21.0", 
                                  "statsmodels",
                                  "git+https://github.com/sjdavenport/pyperm.git"),                     
                     ignore_installed = TRUE)
}
py_config()

test_that("Consistency of lm_test with pyperm package", {
  skip_if_not(py_module_available("pyperm"), "Module pyperm not available")
  
  use_virtualenv("env_test", required = TRUE)
  pyperm <- import("pyperm")
  numpy <- import("numpy")
  
  N <- 100
  DD <- c(1, 10)
  PP <- c(2, 3)
  for (D in DD) {
    for (p in PP) {
      
      X <- matrix(0, nrow = N, ncol = p)
      X[,1] <- 1
      X[,-1] <- runif(N*(p - 1), min = 0, max = 3)
      beta <- matrix(0, nrow = p, ncol = D)
      epsilons <- matrix(rnorm(N*D), nrow = N, ncol = D)
      Y <- X %*% beta + epsilons
      C <- matrix(diag(p)[1, ], nrow = 1)
      
      tstats_py <- pyperm$contrast_tstats_noerrorchecking(lat_data = pyperm$make_field(t(Y)), 
                                                          design = X, 
                                                          contrast_matrix = C)
      tstats_R <- lm_test(Y = Y, X = X, C = C)
      
      expect_equal(t(tstats_py[[1]]$field), 
                   tstats_R$stat_test, 
                   tolerance = 1e-14)
      expect_equal(as.numeric(t(tstats_py[[2]])), 
                   as.numeric(tstats_R$epsilon_est), 
                   tolerance = 1e-14)
      
    }
  }
})

test_that('R/python consistency of `bootstrap_permutation`', {
  skip_if_not(py_module_available("pyperm"), "pyperm module not available")
  
  use_virtualenv("env_test", required = TRUE)
  py_config()
  # load python modules
  pyperm <- import("pyperm")
  numpy <- import("numpy")
  
  N = 500
  DD <- c(2,10)
  PP <- c(2,10)
  LL <- c(1,5)
  
  configs <- expand.grid(
    DD =DD,
    PP = PP,
    LL = LL
  )
  configs <- configs[configs$LL <= configs$PP,]
  seq_configs <- seq_len(nrow(configs))
  
  for (cc in seq_configs) {
    config <- configs[cc, ]
    D <- config[["DD"]]
    p <- config[["PP"]]
    L <- config[["LL"]]
    
    
    X <- matrix(0,nrow = N, ncol = p)
    X[,1] <- 1
    X[, -1] <- runif(N*(p - 1), min = 0, max = 3)
    beta <- matrix(0, nrow = p, ncol = D)
    # Sigma = diag(D) # Need correlation 
    # Sigma = matrix(1, nrow = D, ncol = D) #too much correlation but works for tolerance of 1e-10
    Sigma = toeplitz(0.9^(0:(D-1))) # correlation (similar to real data)
    Sigma_tri <- chol(Sigma)
    Z <- matrix(rnorm(N * D), nrow = N, ncol = D)
    epsilons <- Z %*% t(Sigma_tri)
    Y <- X %*% beta + epsilons
    C <- matrix(diag(p)[1:L, ], nrow = L)
    B <- 500L
    
    res_py <- pyperm$boot_contrasts(lat_data = t(Y), design = X, 
                                    contrast_matrix = C,
                                    n_bootstraps = B, store_boots = TRUE)
    pval_perm <- bootstrap_permutation(Y = Y, X = X, C = C, B, 
                                       replace = TRUE, 
                                       alternative = "two.sided" )
    pval_perm_R <- matrix(pval_perm, nrow = D*nrow(C), ncol = B)
    
    pval_perm_py <- t(res_py[[4]])
    
    m <- D*L
    sample_m <- sample(1:m, min(m, 20), replace = FALSE)
    for (mm in sample_m) {
      expect_true(ks.test(pval_perm_R[mm,], "punif")$p.value > 1e-4)
      expect_true(ks.test(pval_perm_py[mm,], pval_perm_R[mm,])$p.value > 1e-4)
    }
    
    for (mm in sample_m) {
      for (mm2 in sample_m[which(sample_m != mm)]) {
        expect_equal(cor(pval_perm_py[mm,], pval_perm_py[mm2,]),
                     cor(pval_perm_R[mm,], pval_perm_R[mm2,]),
                     tolerance = 0.25)
      }
    }
    
  }
  
  ## Using real dataset 
  
  data(expr_ALL, package = "sanssouci.data")
  groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
  D <- 100
  Y <- expr_ALL[sample(1:nrow(expr_ALL), D), ]
  Y <- expr_ALL[order(matrixStats::rowVars(expr_ALL), decreasing = TRUE)[1:D], ]
  X <- matrix(1,nrow = dim(expr_ALL)[2], ncol = 2)
  X[,2] <- groups
  C <- diag(2)
  B <- 500L
  
  res_py <- pyperm$boot_contrasts(lat_data = Y, design = X, contrast_matrix = C,
                                  n_bootstraps = B, store_boots = TRUE)
  pval_perm <- bootstrap_permutation(Y = t(Y), X = X, C = C, B, 
                                     replace = TRUE, 
                                     alternative = "two.sided" )
  pval_perm_R <- matrix(pval_perm, nrow = nrow(Y)*nrow(C), ncol = B)
  
  pval_perm_py <- t(res_py[[4]])
  
  m <- nrow(Y)*nrow(C)
  sample_m <- sample(1:m, min(m, 50), replace = FALSE)
  for (mm in sample_m) {
    expect_true(ks.test(pval_perm_R[mm,], "punif")$p.value > 1e-4)
    expect_true(ks.test(pval_perm_py[mm,], pval_perm_R[mm,])$p.value > 1e-4)
  }
  
  df <- data.frame()
  for (mm in sample_m) {
    for (mm2 in sample_m[which(sample_m != mm)]) {
      expect_equal(cor(pval_perm_py[mm,], pval_perm_py[mm2,]),
                   cor(pval_perm_R[mm,], pval_perm_R[mm2,]),
                   tolerance = 0.25)
    }
  }
})

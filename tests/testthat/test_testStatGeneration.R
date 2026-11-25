context("Generation of test statistics -- low-level functions")

test_that("Simulate equi-correlated null hypotheses", {
  m <- 1331
  n <- 111
  ## equi-correlated
  rho <- 0.2
  Y <- sim_null_gaussian_equi(m, n, rho)
  expect_equal(nrow(Y), m)
  expect_equal(ncol(Y), n)
})

test_that("Simulation of null hypotheses from factor model", {
  m <- 10
  ## independent
  withr::with_preserve_seed({
    sim0 <- sim_null_gaussian_factor_sv(m)
  })
  expect_equal(length(sim0$Y), m)
  expect_null(sim0$w)
  S0 <- getFactorModelCovarianceMatrix(m)
  expect_identical(S0, diag(m))

  Y0 <- sim_null_gaussian_factor(m, flavor = "independent")
  expect_identical(Y0[, 1], sim0$Y)

  # expected error: "t(P) %*% P should be orthonormal"
  # (but encoding issue => can't get the exact string)
  expect_error(
    sim_null_gaussian_factor_sv(m, c(0, 1))
  )

  ## equi-correlated
  rho <- 0.2
  h <- 1
  P <- matrix(1, m, length(h))
  withr::with_preserve_seed({
    sim1 <- sim_null_gaussian_factor_sv(m, h, P,
      rho = rho
    )
  })
  expect_equal(length(sim1$Y), m)
  expect_equal(length(sim1$W), 1)
  Y1 <- sim_null_gaussian_factor(m,
    flavor = "equi-correlated",
    rho = rho, cov = TRUE
  )
  expect_identical(Y1[, 1], sim1$Y)

  ## 3-factor model
  m <- 4 * floor(m / 4) ## make sure m/4 is an integer
  rho <- 0.5
  h <- c(0.5, 0.3, 0.2) * m
  gamma1 <- rep(1, m)
  gamma2 <- rep(c(-1, 1), each = m / 2)
  gamma3 <- rep(c(-1, 1), times = m / 2)
  P <- cbind(gamma1, gamma2, gamma3) / sqrt(m)

  withr::with_preserve_seed({
    sim3 <- sim_null_gaussian_factor_sv(m, h, P,
      rho = rho
    )
  })
  expect_equal(length(sim3$Y), m)
  expect_equal(length(sim3$W), length(h))

  Y3 <- sim_null_gaussian_factor(m,
    flavor = "3-factor",
    rho = rho
  )
  expect_identical(Y3[, 1], sim3$Y)

  Sigma <- getFactorModelCovarianceMatrix(m, h, P = P, rho = rho)
  expect_equal(nrow(Sigma), ncol(Sigma))
  expect_equal(nrow(Sigma), m)
})

test_that("Simulation of Gaussian null hypotheses from covariance matrix", {
  library("Matrix")
  m <- 100
  n <- 1000

  ## Toeplitz, short range
  tcoefs <- toeplitz((1:m)^(-2))
  Sigma <- Matrix(tcoefs, sparse = TRUE)
  Y <- sim_null_gaussian_sigma(Sigma, n)
  SigmaHat <- Y %*% t(Y) / n
  expect_equal(nrow(Y), m)
  expect_equal(ncol(Y), n)
})

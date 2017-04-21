context("Generation of test statistics")


test_that("Simulate equi-correlated null hypotheses", {
    m <- 1331
    n <- 111
    ## equi-correlated
    rho <- 0.2
    Y <- simulateGaussianEquiCorrelatedNulls(m, n, rho)
    expect_equal(nrow(Y), m)
    expect_equal(ncol(Y), n)
})

test_that("Simulate Gaussian equi-correlated test statistics", {
    m <- 123
    rho <- 0.2
    B <- 100
    pi0 <- 0.5
    sim <- simulateEqui(m, rho, B, pi0, SNR=1)
    expect_equal(length(sim$x), m)
    expect_equal(nrow(sim$X0), m)
    expect_equal(length(sim$H), m)
    expect_equal(ncol(sim$X0), B)
    expect_equal(round(m*pi0), m-sum(sim$H))
    
    stat <- sim$x
    H <- sim$H
    sa <- tapply(stat, H, mean)
    expect_gt(sa[2], sa[1])  ## test stat is larger under H1
    
})


test_that("Simulation of null hypotheses from factor model", {
    m <- 10
    ## independent
    sim0 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m)
    expect_equal(length(sim0$Y), m)
    S0 <- getFactorModelCovarianceMatrix(m)
    expect_identical(S0, diag(m))

    ## equi-correlated
    rho <- 0.2
    h <- 1
    P <- matrix(1, m, length(h))
    sim1 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
    expect_equal(length(sim1$Y), m)
    
    ## 3-factor model
    m <- 4*floor(m/4) ## make sure m/4 is an integer
    rho <- 0.5
    h <- c(0.5, 0.3, 0.2)*m
    gamma1 <- rep(1,m)
    gamma2 <- rep(c(-1, 1), each=m/2)
    gamma3 <- rep(c(-1, 1), times=m/2)
    P <- cbind(gamma1, gamma2, gamma3)/sqrt(m)
    
    sim3 <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)
    expect_equal(length(sim3$Y), m)
    expect_equal(length(sim3$W), length(h))
})

test_that("Simulation of Gaussian null hypotheses from covariance matrix", {
    
    library("Matrix")
    m <- 100
    n <- 1000

    ## Toeplitz, short range
    tcoefs <- toeplitz((1:m)^(-2))
    Sigma <- Matrix(tcoefs, sparse = TRUE)
    Y <- simulateGaussianNullsFromSigma(Sigma, n)
    SigmaHat <- Y %*% t(Y)/n
    expect_equal(nrow(Y), m)
    expect_equal(ncol(Y), n)
})

test_that("Simulate test statistics as in Meinshausen (2006)", {
    m <- 123
    rho <- 0.2
    n <- 100
    pi0 <- 0.5
    B <- 1e3

    sim <- simulateMein2006(m, rho, n, B, pi0, SNR=1)
    expect_equal(length(sim$x), m)
    expect_equal(nrow(sim$X0), m)
    expect_equal(length(sim$H), m)
    expect_equal(ncol(sim$X0), B)
    expect_equal(round(m*pi0), m-sum(sim$H))

    stat <- sim$x
    H <- sim$H
    sa <- tapply(stat, H, mean)
    expect_gt(sa[2], sa[1])  ## test stat is larger under H1
})


test_that("Simulate Gaussian test statistics with Toeplitz covariance", {
    m <- 123
    pow <- -0.5
    B <- 100
    pi0 <- 0.5

    sim <- simulateToeplitz(m, pow, B, pi0, SNR=1)
    stat <- sim$x
    H <- sim$H
    sa <- tapply(stat, H, mean)
    expect_gt(sa[2], sa[1])  ## test stat is larger under H1
})

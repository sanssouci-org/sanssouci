context("Simulation of Gaussian samples")

m <- 123
rho <- 0.2
n <- 100
pi0 <- 0.5

test_that("output dimensions are as expected for two-sample tests", {
    sim <- gaussianSamples(m, rho, n, pi0, SNR = 1, prob = 0.5)
    
    expect_equal(ncol(sim$X), n)
    expect_equal(nrow(sim$X), m)
    expect_length(sim$H, m)
    
    tbl <- table(colnames(sim$X))
    expect_length(tbl, 2)
    
    tbl <- table(sim$H)
    expect_length(tbl, 2)
})

test_that("output dimensions are as expected for one-sample tests", {
    sim <- gaussianSamples(m, rho, n, pi0, SNR = 1, prob = 1)
    
    expect_equal(ncol(sim$X), n)
    expect_equal(nrow(sim$X), m)
    expect_length(sim$H, m)
    
    tbl <- table(colnames(sim$X))
    expect_length(tbl, 1)
    
    tbl <- table(sim$H)
    expect_length(tbl, 2)
})
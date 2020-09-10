context("Simulation of Gaussian test statistics")

test_that("output dimensions are as expected", {
    m <- 123
    B <- 100
    sim <- gaussianTestStatistics(m, B)
    
    expect_length(sim$x, m)
    expect_equal(nrow(sim$X0), m)
    expect_equal(ncol(sim$X0), B)
    expect_length(sim$H, m)
  
    pi0 <- 0.8
    m0 <- round(m * pi0)
    test <- replicate(10, {
        sim <- gaussianTestStatistics(m, B, pi0 = pi0)
        expect_equal(m0, m - sum(sim$H))
    })
})

test_that("catching parameter incompatibilities", {
    m <- 123
    B <- 100
    expect_error(gaussianTestStatistics(m, B, dep = "equi", param = -1))
    expect_error(gaussianTestStatistics(m, B, param = -1))
    expect_error(gaussianTestStatistics(m, B, dep = "Toeplitz", param = 1))
})

test_that("mean under the alternative larger than under the null", {
    m <- 123
    B <- 100
    pi0 <- 0.8

    sim <- gaussianTestStatistics(m, B, pi0 = pi0, dep = "equi", param = 0, SNR = 1)
    sa <- tapply(sim$x, sim$H, mean)
    expect_gt(sa[2], sa[1])
    
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, dep = "equi", param = 0.3, SNR = 1)
    sa <- tapply(sim$x, sim$H, mean)
    expect_gt(sa[2], sa[1])  
    
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, dep = "Toeplitz", param = -0.5, SNR = 1)
    sa <- tapply(sim$x, sim$H, mean)
    expect_gt(sa[2], sa[1]) 
})



test_that("Argument 'SNR'", {
    m <- 123
    B <- 100
    pi0 <- 0.8
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = 1)

    m1 <- m - round(m * pi0)
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = rep(1, m1))

    expect_error(gaussianTestStatistics(m, B, pi0 = pi0, SNR = rep(1, 2*m)))
})

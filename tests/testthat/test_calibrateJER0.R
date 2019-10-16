context("Calibration of reference family")

test_that("Numerical reproducibility of one particular call to 'calibrateJER0'", {
    set.seed(0xBEEF)
    m <- 3210
    B <- 123
    rho <- 0
    pi0 <- 0.9
    SNR <- 2
    
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = SNR, dep = "equi", param = rho)
    X0 <- sim$X0
    x <- sim$x
    
    res <- calibrateJER0(X0, refFamily = "kFWER", alpha = 0.2, stat = x, kMax = m, verbose = FALSE)
    
    pathname <- system.file("extdata/calibrateJER0_results.rds", package = "sansSouci")
    ref <- readRDS(pathname)
    
    expect_equal(ref, res)
})

context("Calibration of reference family")

test_that("Numerical reproducibility of one particular call to 'calibrateJER0'", {
    set.seed(0xBEEF)
    m <- 3210
    B <- 123
    rho <- 0
    pi0 <- 0.9
    SNR <- 3
    
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = SNR, dep = "equi", param = rho)
    pval <- pnorm(sim$x, lower.tail = FALSE)
    p0 <- pnorm(sim$X0, lower.tail = FALSE)
    res <- calibrateJER0(mat = p0, refFamily = "kFWER", alpha = 0.5, p.values = pval, kMax = m, verbose = FALSE)

    pathname <- system.file("extdata/calibrateJER0_results.rds", package = "sansSouci")
    ref <- readRDS(pathname)
    
    expect_equal(ref, res)
})

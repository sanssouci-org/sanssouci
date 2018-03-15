context("Estimation of a joint FWER controlling family by lambda-adjustment")

test_that("Numerical reproducibility of one particular call to 'jointFWERControl'", {
    set.seed(0xBEEF)
    m <- 1000
    B <- 1000
    rho <- 0
    pi0 <- 0.9
    SNR <- 2
    
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = SNR, dep = "equi", param = rho)
    X0 <- sim$X0
    x <- sim$x
    
    res <- jointFWERControl(X0, refFamily="kFWER", alpha=0.2, stat=x, kMax=m, verbose=FALSE)
    
    pathname <- system.file("extdata/jointFWERControlResults.rds", package="sansSouci")
    ref <- readRDS(pathname)
    
    expect_equal(ref, res)
})

context("JER calibration")

test_that("JER calibration -- vanilla tests", {
    m <- 123
    rho <- 0.2
    n <- 100
    pi0 <- 0.5
    sim <- gaussianSamples(m, rho, n, pi0, SNR = 2, prob = 0.5)
    X <- sim$X
    alpha <- 0.2
    B <- 100

    cal <- calibrateJER(X, B, alpha, refFamily = "Simes")
    expect_length(cal$thr, m)
    expect_length(cal$stat, m)

    cal <- calibrateJER(X, B, alpha, refFamily = "Simes", K = 50)
    expect_length(cal$thr, 50)
    expect_length(cal$stat, m)
})



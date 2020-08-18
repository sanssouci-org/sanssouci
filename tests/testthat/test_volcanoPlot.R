test_that("Vanilla test for 'volcanoPlot'", {
    m <- 500
    pi0 <- 0.5
    m1 <- m-m*pi0
    sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
                           pi0 = pi0, SNR = runif(m1)*6-3, prob = 0.5)
    X <- sim$X
    alpha <- 0.2
    cal <- calibrateJER(X, B = 1e2, alpha = alpha, refFamily="Simes")
    vp <- volcanoPlot(X, categ = colnames(X), thr = cal$thr, p = 3, r = 0.3, ylim = c(0, 6))
    expect_null(vp)
})
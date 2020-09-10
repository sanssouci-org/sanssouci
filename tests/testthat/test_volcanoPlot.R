context("Volcano plots")

test_that("Vanilla test for 'volcanoPlot'", {
    m <- 500
    pi0 <- 0.5
    m1 <- m-m*pi0
    SNR <- 5*(runif(m1)-0.5)
    sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
                           pi0 = pi0, SNR = SNR, prob = 0.5)
    X <- sim$X
    categ <- sim$categ
    alpha <- 0.2
    cal <- calibrateJER(X, categ, B = 1e2, alpha = alpha, refFamily="Simes")
    
    vp <- volcanoPlot(X, categ, thr = cal$thr, 
                      p = 1, q = 1, r = 0)
    expect_equal(vp, seq_len(m))     ## no active filter: all genes selected
    
    vp <- volcanoPlot(X, categ, thr = cal$thr, 
                      p = 1, q = 1, r = Inf)
    expect_equal(length(vp), 0)      ## too stringent filtering on fold change

    vp <- volcanoPlot(X, categ, thr = cal$thr, 
                      p = 0, q = 1, r = 0)
    expect_equal(length(vp), 0)      ## too stringent filtering on p-value
    
    vp <- volcanoPlot(X, categ, thr = cal$thr, 
                      p = 1, q = 0, r = 0)
    expect_equal(length(vp), 0)      ## too stringent filtering on q-value

    ## filtering both on p-values and q-values
    expect_warning(volcanoPlot(X, categ, thr = cal$thr, 
                               p = 0.01, q = 0.05, r = 0, ylim = c(0, 6)))
})
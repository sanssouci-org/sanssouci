context("post-hoc bounds derived from Simes' inequality")

require_cherry <- function() {
    if (!require("cherry", quietly = TRUE)) {
        skip("Package 'cherry' not available")
    }
}

test_cherry <- function() {
    skip("Test not working yet!")
}

test_that("posthocBySimes is cherry::pickSimes on 'NAEP' data", {
    require_cherry()
    # test_cherry()

    data("NAEP", package="cherry")
    m <- length(NAEP)
    nms <- names(NAEP)
    hom <- cherry::hommelFast(NAEP)
    p <- NAEP
    
    test <- replicate(10, {
        nR <- round(runif(1)*m)
        R <- sample(nms, nR)
        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)
        res_BNR <- posthocBySimes(p, R, alpha = alpha)
        expect_equal(res_cherry, res_BNR) 
        
        # artisanal way
        pi0_hat <-  1 - posthocBySimes0(p, 1:m, alpha = alpha)/m  ## in theory one should iterate this
        res_BNR0 <- posthocBySimes0(p, R, alpha = alpha/pi0_hat)
        expect_equal(res_cherry, res_BNR0) 
    })
})

test_that("posthocBySimes is cherry::pickSimes on simulated data", {
    require_cherry()
    # test_cherry()
    
    m <- 1e2
    m1 <- 10
    
    test <- replicate(10, {
        p <- 1 - pnorm(c(rnorm(m1, mean = 4), rnorm(m - m1, mean = 0)))
        hom <- cherry::hommelFast(p)
        nR <- round(runif(1)*m)
        R <- sample(nR)

        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)
        res_BNR <- posthocBySimes(p, R, alpha = alpha)
        expect_equal(res_cherry, res_BNR) 
        
        # artisanal way
        pi0_hat <-  1 - posthocBySimes0(p, 1:m, alpha = alpha)/m  ## in theory one should iterate this
        res_BNR0 <- posthocBySimes0(p, R, alpha = alpha/pi0_hat)
        expect_equal(res_cherry, res_BNR0) 
    })
})

test_that("posthocBySimes is cherry::pickSimes on larger simulated data", {
    require_cherry()
    # test_cherry()
    
    m <- 1e3
    m1 <- 100
    
    test <- replicate(10, {
        p <- 1 - pnorm(c(rnorm(m1, mean = 4), rnorm(m - m1, mean = 0)))
        hom <- hommelFast(p)
        nR <- round(runif(1)*m)
        R <- sample(nR)

        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)
        res_BNR <- posthocBySimes(p, R, alpha = alpha)
        expect_equal(res_cherry, res_BNR) 
        
        # artisanal way
        pi0_hat <-  1 - posthocBySimes0(p, 1:m, alpha = alpha)/m  ## in theory one should iterate this
        res_BNR0 <- posthocBySimes0(p, R, alpha = alpha/pi0_hat)
        expect_equal(res_cherry, res_BNR0) 
    })
})


test_that("posthocBySimes0 is posthocBySimes0Rcpp for simulated data", {
    m <- 1e3
    m1 <- 100
    alphas <- seq(from = 0, to = 1, by = 0.1)
    p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
    nR <- round(runif(1)*m)
    R <- sample(nR)
    for (alpha in alphas) {
        expect_equal(posthocBySimes0(p, R, alpha = alpha),
                     posthocBySimes0Rcpp(p, R, alpha = alpha))
    }
})

test_that("posthocBySimes is posthocBySimes0Rcpp for simulated data an no signal", {
    m <- 1e3
    m1 <- 100
    alpha <- 0.01
    p <- 1-pnorm(c(rnorm(m1, mean=0), rnorm(m-m1, mean=0)))
    nR <- round(runif(1)*m)
    R <- sample(nR)
    expect_equal(posthocBySimes(p, R, alpha = alpha),
                 posthocBySimes0Rcpp(p, R, alpha = alpha))
    expect_equal(posthocBySimes(p, R, alpha = alpha, Rcpp = TRUE),
                 posthocBySimes0Rcpp(p, R, alpha = alpha))
})

test_that("posthocBySimes0 can be reproduced by minTP", {
    m <- 1e3
    m1 <- 100
    alphas <- seq(from = 0, to = 1, by = 0.1)
    
    x <- c(rnorm(m1, mean=4), rnorm(m-m1, mean=0))
    sx <- sort(x, decreasing = TRUE)
    p <- 1 - pnorm(x)
    nR <- round(runif(1)*m)
    R <- sample(nR)
    for (alpha in alphas) {
        thrSimes <- SimesThresholdFamily(m)(alpha)
        ubSimes <- minTP(p[R], thrSimes)
        expect_equal(posthocBySimes0(p, R, alpha = alpha), 
                     ubSimes)
    }
    
    posthocBySimes0(p, R, alpha = alpha)
})
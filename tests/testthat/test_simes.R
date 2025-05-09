context("post-hoc bounds derived from Simes' inequality")

require_cherry <- function() {
    if (!require("cherry", quietly = TRUE)) {
        skip("Package 'cherry' not available")
    }
}

test_cherry <- function() {
    skip("Test not working yet!")
}

test_that("BNR+Simes is cherry::pickSimes on 'NAEP' data", {
    require_cherry()
    # test_cherry()

    data("NAEP", package="cherry")
    m <- length(NAEP)
    nms <- names(NAEP)
    hom <- cherry::hommelFast(NAEP)
    p <- NAEP
    
    test <- replicate(10, {
        nR <- ceiling(runif(1)*m)
        R <- sample(nms, nR)
        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)

        thr <- t_linear(alpha, 1:m, m)
        m0_hat <- maxFP(p.values = p, thr = thr)
        thr <- t_linear(alpha, 1:m0_hat, m0_hat)
        res_BNR <- minTP(p[R], thr = thr)
        expect_equal(res_cherry, res_BNR)
    })
})

test_that("BNR+Simes is cherry::pickSimes on simulated data", {
    require_cherry()
    # test_cherry()
    
    m <- 1e2
    m1 <- 10
    
    test <- replicate(10, {
        p <- 1 - pnorm(c(rnorm(m1, mean = 3), rnorm(m - m1, mean = 0)))
        hom <- cherry::hommelFast(p)
        nR <- ceiling(runif(1)*m)
        R <- sample(nR)

        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)

        thr <- t_linear(alpha, 1:m, m)
        m0_hat <- maxFP(p.values = p, thr = thr)
        thr <- t_linear(alpha, 1:m0_hat, m0_hat)
        res_BNR <- minTP(p[R], thr = thr)
        expect_equal(res_cherry, res_BNR)
    })
})

test_that("BNR+Simes is cherry::pickSimes on larger simulated data", {
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

        thr <- t_linear(alpha, 1:m, m)
        m0_hat <- maxFP(p.values = p, thr = thr)
        thr <- t_linear(alpha, 1:m0_hat, m0_hat)
        res_BNR <- minTP(p[R], thr = thr)
        expect_equal(res_cherry, res_BNR)
    })
})
# 
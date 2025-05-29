context("post-hoc bounds derived from Simes' inequality")

require_cherry <- function() {
    if (!require("cherry", quietly = TRUE)) {
        skip("Package 'cherry' not available")
    }
}

test_that("BNR+Simes is cherry::pickSimes on 'NAEP' data", {
    require_cherry()

    data("NAEP", package="cherry")
    m <- length(NAEP)
    hom <- cherry::hommelFast(NAEP)
    p <- NAEP
    
    test <- replicate(10, {
        nR <- ceiling(runif(1)*m)
        R <- sample(nR)
        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)
        res_BNR <- posthocBySimes(p, R, alpha, flavor = "full step down")
        expect_equal(res_cherry, res_BNR)
    })
})

test_that("BNR+Simes is cherry::pickSimes on simulated data", {
    require_cherry()

    m <- 1e2
    m1 <- 10
    
    test <- replicate(10, {
        p <- 1 - pnorm(c(rnorm(m1, mean = 3), rnorm(m - m1, mean = 0)))
        hom <- cherry::hommelFast(p)
        nR <- ceiling(runif(1)*m)
        R <- sample(nR)

        alpha <- 0.1
        res_cherry <- cherry::pickSimes(hom, R, alpha = alpha, silent=TRUE)
        res_BNR <- posthocBySimes(p, R, alpha)
        expect_equal(res_cherry, res_BNR)
    })
})

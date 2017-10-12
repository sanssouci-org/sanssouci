context("post-hoc bounds derived from Simes' inequality")

require_cherry <- function() {
    if (!require("cherry")) {
        skip("Package 'cherry' not available")
    }
}

test_that("posthocBySimes is cherry::pickSimes on 'NAEP' data", {
    require_cherry()

    data("NAEP", package="cherry")
    m <- length(NAEP)
    nms <- names(NAEP)
    hom <- cherry::hommelFast(NAEP)
    
    test <- replicate(10, {
        nR <- round(runif(1)*m)
        R <- sample(nms, nR)
        expect_equal(cherry::pickSimes(hom, R, silent=TRUE),
                     posthocBySimes(NAEP, R))
    })
})

test_that("posthocBySimes is cherry::pickSimes on simulated data", {
    require_cherry()
    
    m <- 1e2
    m1 <- 10
    
    test <- replicate(10, {
        p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
        hom <- cherry::hommelFast(p)
        nR <- round(runif(1)*m)
        R <- sample(nR)
        expect_equal(cherry::pickSimes(hom, R, silent=TRUE),
                     posthocBySimes(p, R))
    })
})

test_that("posthocBySimes is cherry::pickSimes on larger simulated data", {
    require_cherry()
    m <- 1e3
    m1 <- 100
    
    test <- replicate(10, {
        p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
        hom <- hommelFast(p)
        nR <- round(runif(1)*m)
        R <- sample(nR)
        expect_equal(pickSimes(hom, R, silent=TRUE),
                     posthocBySimes(p, R), tolerance=1)
    })
})


test_that("posthocBySimes is posthocBySimesRcpp for simulated data", {
    m <- 1e3
    m1 <- 100
    alphas <- seq(from = 0, to = 1, by = 0.1)
    test <- replicate(10, {
        p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
        nR <- round(runif(1)*m)
        R <- sample(nR)
        for (alpha in alphas) {
            expect_equal(posthocBySimes(p, R, alpha = alpha),
                         posthocBySimesRcpp(p, R, alpha = alpha))
        }
    })
})

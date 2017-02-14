library("sansSouci")


context("post-hoc bounds derived from Simes' inequality")

test_that("posthocBySimes is cherry::pickSimes on 'NAEP' data", {
    if (require("cherry")) {
        ## Data from package `cherry`
        data("NAEP", package="cherry")
        m <- length(NAEP)
        nms <- names(NAEP)
        hom <- hommelFast(NAEP)
        
        ntests <- 10
        for (tt in 1:ntests) {
            nR <- round(runif(1)*m)
            R <- sample(nms, nR)
            expect_equal(pickSimes(hom, R, silent=TRUE),
                         posthocBySimes(NAEP, R))
        }
    }
})

test_that("posthocBySimes is cherry::pickSimes on simulated data", {
    skip("Not working yet!")
    if (require("cherry")) {
        m <- 1e2
        m1 <- 10
        
        ntests <- 10
        for (tt in 1:ntests) {
            p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
            hom <- hommelFast(p)
            nR <- round(runif(1)*m)
            R <- sample(nR)
            expect_equal(pickSimes(hom, R, silent=TRUE),
                         posthocBySimes(p, R))
        }
    }
})

test_that("posthocBySimes is cherry::pickSimes on larger simulated data", {
    if (require("cherry")) {
        m <- 1e3
        m1 <- 100
        
        ntests <- 10
        for (tt in 1:ntests) {
            p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
            hom <- hommelFast(p)
            nR <- round(runif(1)*m)
            R <- sample(nR)
            expect_equal(pickSimes(hom, R, silent=TRUE),
                         posthocBySimes(p, R), tolerance=1)
        }
    }
})


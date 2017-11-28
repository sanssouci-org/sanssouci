context("Calculation of permutation p-values")

p <- 513
n <- 38
mat <- matrix(rnorm(p*n), ncol = n)


cls <- rep(c(0, 1), times = c(27, n - 27))
res <- testByTwoSamplePermutation(X = mat, cls = cls, B = 1000)

test_that("Two different ways of getting permutation p-values give identical results", {
    
    T0 <- cbind(res$T0, res$T)
    sw <- sweep(abs(T0), MARGIN=1, STATS=abs(res$T), FUN=">=")
    p <- rowMeans(sw)
    expect_identical(res$p, p)
})

test_that("Correctness of permutation p-values", {

    cls <- rep(c(0, 1), times = c(27, n-27))
    res <- testByTwoSamplePermutation(X = mat, cls = cls, B = 1000)
    
    ## stochastic domination by Uniform distribution
    for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
        cs <- colSums(res$p0<=alpha)
        expect_true(mean(cs/p)<=alpha)  
    }    
})

res <- testBySignFlipping(X = mat, B = 1000)

test_that("Two different ways of getting sign flipping p-values give identical results", {
    T0 <- cbind(res$T0, res$T)
    sw <- sweep(abs(T0), MARGIN=1, STATS=abs(res$T), FUN=">=")
    p <- rowMeans(sw)
    expect_identical(res$p, p)
})


test_that("Correctness of sign-flipping p-values", {
    
    ## stochastic domination by Uniform distribution
    for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
        cs <- colSums(res$p0<=alpha)
        expect_true(mean(cs/p)<=alpha)  
    }    
})

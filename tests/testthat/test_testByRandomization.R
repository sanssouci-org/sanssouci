context("Calculation of randomization p-values")

p <- 513
n <- 38
mat <- matrix(rnorm(p*n), ncol = n)


cls <- rep(c(0, 1), times = c(27, n - 27))
res <- testByRandomization(X = mat, B = 1000, cls = cls)
expect_equal(res$flavor, "perm")
rm(cls)

test_that("Two different ways of getting permutation p-values give identical results", {
    
    T0 <- cbind(res$T0, res$T)
    sw <- sweep(abs(T0), MARGIN=1, STATS=abs(res$T), FUN=">=")
    p <- rowMeans(sw)
    expect_identical(res$p, p)
})

test_that("Correctness of permutation p-values", {

    ## stochastic domination by Uniform distribution
    for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
        cs <- colSums(res$p0<=alpha)
        expect_true(mean(cs/p)<=alpha)  
    }    
})

res <- testByRandomization(X = mat, B = 1000)
expect_equal(res$flavor, "flip")

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

test_that("Consistency of Rcpp and R sign-flipping p-values", {
    m <- 123
    n <- 38
    B <- 100
    
    X <- matrix(rnorm(m*n), ncol=n)

    set.seed(123)
    T <- testBySignFlipping(X, B)
    
    set.seed(123)
    TR <- testBySignFlippingR(X, B)
    expect_equal(T, TR)
})



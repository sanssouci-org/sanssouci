context("Calculation of randomization p-values")

m <- 513
n <- 38
mat <- matrix(rnorm(m*n), ncol = n)

test_that("Correctness of permutation p-values", {
    cls <- rep(c(0, 1), times = c(27, n - 27))
    res <- testByRandomization(X = mat, B = 1000, cls = cls, rand.p.value = TRUE)
    expect_equal(res$flavor, "perm")
    rm(cls)

    # consistency between parametric and randomization p-values
    expect_gt(cor(res$p, res$rand.p), 0.9)
    
    # Two different ways of getting permutation p-values
    T0 <- cbind(res$T0, res$T)
    sw <- sweep(abs(T0), MARGIN = 1, STATS = abs(res$T), FUN = ">=")
    p <- rowMeans(sw)
    expect_identical(res$rand.p, p)

    ## stochastic domination by Uniform distribution
    for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
        cs <- colSums(res$rand.p0 <= alpha)
        expect_true(mean(cs/m) <= alpha)  
    }    
})


test_that("Correctness of sign-flipping p-values", {
    res <- testByRandomization(X = mat, B = 1000, rand.p.value = TRUE)
    expect_equal(res$flavor, "flip")

    # Two different ways of getting sign flipping p-values give identical results"
    T0 <- cbind(res$T0, res$T)
    sw <- sweep(abs(T0), MARGIN = 1, STATS = abs(res$T), FUN = ">=")
    p <- rowMeans(sw)
    expect_identical(res$rand.p, p)

    # consistency between parametric and randomization p-values
    for (alt in c("two.sided", "less", "greater")) {
        res <- testByRandomization(X = mat, B = 1000, rand.p.value = TRUE, alternative = alt)
        expect_gt(cor(res$p, res$rand.p), 0.9)
    }
    
    ## stochastic domination by Uniform distribution
    for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
        cs <- colSums(res$rand.p0 <= alpha)
        expect_true(mean(cs/m) <= alpha)  
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


test_that("Consistency of randomization p-values with different alternatives", {
    for (alt in c("two.sided", "less", "greater")) {
        
        ## permutation
        cls <- rep(c(0, 1), times = c(27, n - 27))
        two.sided <- testByRandomization(X = mat, B = 100, cls = cls, alternative = "two.sided")
        greater <- testByRandomization(X = mat, B = 100, cls = cls, alternative = "greater")
        less <- testByRandomization(X = mat, B = 100, cls = cls, alternative = "less")
        expect_equal(two.sided$p, 2*pmin(greater$p, less$p))

        ## swapping the class labels        
        clsInv <- 1 - cls
        two.sidedInv <- testByRandomization(X = mat, B = 100, cls = clsInv, alternative = "two.sided")
        greaterInv <- testByRandomization(X = mat, B = 100, cls = clsInv, alternative = "greater")
        lessInv <- testByRandomization(X = mat, B = 100, cls = clsInv, alternative = "less")
        expect_equal(two.sidedInv$p, two.sided$p)
        expect_equal(greaterInv$p, less$p)
        expect_equal(lessInv$p, greater$p)
        
        ## randomization
        two.sided <- testByRandomization(X = mat, B = 100, alternative = "two.sided")
        greater <- testByRandomization(X = mat, B = 100, alternative = "greater")
        less <- testByRandomization(X = mat, B = 100, alternative = "less")
        expect_equal(two.sided$p, 2*pmin(greater$p, less$p))
    }
})
context("Calculation of Welch test statistics and p-values")

p <- 2e2
n <- 100
## null distribution
mat <- matrix(rnorm(p*n), ncol=n)
cls <- rep(c(0, 1), times=c(n/2, n-n/2))

test_that("rowWelchTests <=> welch.test", {
    for (alt in c("two.sided", "greater", "less")) {
    fwt <- rowWelchTests(mat, categ = cls, alternative = alt)
    
        # Ordinary Welch t-test (reasonably fast for 200 tests)
        wt <- apply(mat, 1, FUN = function(x) {
            tt <- t.test(x[cls == 1], x[cls == 0], alternative = alt)  ## test stat positive if "1 > 0"
            c(tt$p.value, tt$statistic, tt$parameter)
        })
        dwt <- as.data.frame(t(wt))
        names(dwt) <- c("p.value", "statistic", "parameter")
        expect_equal(fwt$p.value, dwt$p.value, tolerance = 1e-12)
        expect_equal(fwt$statistic, dwt$statistic, tolerance = 1e-12)
        expect_equal(fwt$parameter, dwt$parameter, tolerance = 1e-12)
    }
})

test_that("correctness of rowWelchTests alternative", {
    
    p <- 30
    n <- 1000
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol=n)
    cls <- rep(c(0, 1), times=c(n/2, n-n/2))
    
    ## adding some signal
    i1 <- which(cls==1)
    idx_greater <- 1:10
    idx_less <- 11:20
    idx_two.sided <- c(idx_greater, idx_less)
    idx_null <- setdiff(1:p, idx_two.sided)
    mat[idx_greater, i1] <- 1 + mat[idx_greater, i1]  
    mat[idx_less, i1] <- -1 + mat[idx_less, i1]  
    
    alt <- "two.sided"    
    pval <- rowWelchTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mt, m0)
    expect_lt(ml, m0)
    expect_lt(mg, m0)
    
    alt <- "greater"    
    pval <- rowWelchTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mg, m0)
    expect_gt(ml, m0)
    
    alt <- "less"    
    pval <- rowWelchTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_gt(mg, m0)
    expect_lt(ml, m0)
})


test_that("suffWelchTest gives identical results to welch.test", {

    mat <-     
    replicate(10, {
        x <- rnorm(1234)
        y <- rnorm(2345)
        target <- t.test(x, y)
        swt <- suffWelchTest(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
        
        expect_equivalent(swt$statistic, target$statistic)
        expect_equal(swt$p.value, target$p.value)
        expect_equivalent(swt$parameter, target$parameter)
    })
})


test_that("Sanity checks of 'categCheck' throw errors when expected to", {
    expect_error(categCheck(c(1,2), 2), 
                            "Expected two categories named '0' and '1'!")
})

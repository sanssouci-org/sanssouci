context("Calculation of Wilcoxon test statistics and p-values")

test_that("rowWilcoxonTests <=> wilcox.test", {
    
    p <- 2e2
    n <- 100
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol = n)
    cls <- rep(c(0, 1), times = c(n/2, n-n/2))
    
    for (alt in c("two.sided", "greater", "less")) {
        fwt <- rowWilcoxonTests(mat = mat, categ = cls, alternative = alt)
        
        # Ordinary Wilcoxon test (reasonably fast for 200 tests)
        wt <- apply(mat, 1, FUN = function(x) {
            wt <- wilcox.test(x[cls == 1], x[cls == 0], alternative = alt)  ## test stat large if "1 > 0"
            c(wt$p.value, wt$statistic)
        })
        dwt <- as.data.frame(t(wt))
        names(dwt) <- c("p.value", "statistic")
        expect_equal(fwt$p.value, dwt$p.value, tolerance = 1e-12)
        expect_equal(fwt$statistic, dwt$statistic, tolerance = 1e-12)
    }
})


test_that("correctness of rowWilcoxonTests alternative", {
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
    pval <- rowWilcoxonTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mt, m0)
    expect_lt(ml, m0)
    expect_lt(mg, m0)

    alt <- "greater"    
    pval <- rowWilcoxonTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mg, m0)
    expect_gt(ml, m0)
    
    alt <- "less"    
    pval <- rowWilcoxonTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_gt(mg, m0)
    expect_lt(ml, m0)
})


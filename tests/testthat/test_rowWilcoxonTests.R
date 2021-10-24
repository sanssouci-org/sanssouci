context("Calculation of Wilcoxon test statistics and p-values")

test_that("rowWilcoxonTests <=> wilcox.test", {
    p <- 53
    n <- 35
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol = n)
    cls <- rbinom(n, 1, 1/3)

    alt <- sample(c("two.sided", "greater", "less"), 1)
    fwt <- rowWilcoxonTests(mat = mat, categ = cls, alternative = alt)

    # Ordinary Wilcoxon test (reasonably fast for small number of tests)
    wt <- apply(mat, 1, FUN = function(x) {
        wt <- wilcox.test(x[cls == 1], x[cls == 0], alternative = alt, exact = FALSE)  ## test stat large if "1 > 0"
        c(wt$p.value, wt$statistic)
    })
    dwt <- as.data.frame(t(wt))
    names(dwt) <- c("p.value", "statistic")
    expect_equal(fwt$p.value, dwt$p.value)
    expect_equal(fwt$statistic, dwt$statistic)
})


test_that("rowWilcoxonTests <=> wilcox.test (several contrasts at a time)", {
    p <- 15
    n <- 35
    mat_small <- matrix(rnorm(p*n), ncol = n)
    alt <- sample(c("two.sided", "greater", "less"), 1)
    cls_mat <- replicate(10, rbinom(n, 1, 1/3))
    fwt <- rowWilcoxonTests(mat_small, categ = cls_mat, alternative = alt)

    # Ordinary Wilcoxon test (reasonably fast for small number of tests)
    nr <- nrow(mat_small)
    nc <- ncol(cls_mat)
    res <- matrix(NA_real_, nr, nc)
    for (rr in seq_len(nr)) {
        xy <- mat_small[rr, ]
        for (cc in seq_len(nc)) {
            cp <- cls_mat[, cc]
            tt <- wilcox.test(xy[cp == 1], 
                         xy[cp == 0], 
                         alternative = alt, 
                         exact = FALSE)  ## test stat positive if "1 > 0"
            res[rr, cc] <- tt$p.value
        }
    }
    expect_equal(fwt$p.value, res)
})


test_that("Direction of 'alternative' in rowWilcoxonTests", {
#    skip("Test not working yet!")
    p <- 30
    n <- 1000
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol=n)
    cls <- rbinom(n, 1, 1/2)
    
    ## adding some signal
    i1 <- which(cls == 1)
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

### test for rowWicoxonTestsV1
context("Calculation of Wilcoxon test statistics and p-values (V1)")

test_that("rowWilcoxonTestsV1 <=> wilcox.test", {
    p <- 53
    n <- 35
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol = n)
    cls <- rbinom(n, 1, 1/3)
    
    alt <- sample(c("two.sided", "greater", "less"), 1)
    fwt <- rowWilcoxonTestsV1(mat = mat, categ = cls, alternative = alt)
    fwt1 <- rowWilcoxonTests1V1(mat = mat, categ = cls, alternative = alt)
    expect_equal(fwt, fwt1)
    
    # Ordinary Wilcoxon test (reasonably fast for small number of tests)
    wt <- apply(mat, 1, FUN = function(x) {
        wt <- wilcox.test(x[cls == 1], x[cls == 0], alternative = alt, exact = FALSE)  ## test stat large if "1 > 0"
        c(wt$p.value, wt$statistic)
    })
    dwt <- as.data.frame(t(wt))
    names(dwt) <- c("p.value", "statistic")
    expect_equal(fwt$p.value, dwt$p.value)
    expect_equal(fwt$statistic, dwt$statistic)
})


test_that("rowWilcoxonTests <=> wilcox.test (several contrasts at a time)", {
    p <- 15
    n <- 35
    mat_small <- matrix(rnorm(p*n), ncol = n)
    alt <- sample(c("two.sided", "greater", "less"), 1)
    cls_mat <- replicate(10, rbinom(n, 1, 1/3))
    fwt <- rowWilcoxonTestsV1(mat_small, categ = cls_mat, alternative = alt)
    
    # Ordinary Wilcoxon test (reasonably fast for small number of tests)
    nr <- nrow(mat_small)
    nc <- ncol(cls_mat)
    res <- matrix(NA_real_, nr, nc)
    for (rr in seq_len(nr)) {
        xy <- mat_small[rr, ]
        for (cc in seq_len(nc)) {
            cp <- cls_mat[, cc]
            tt <- wilcox.test(xy[cp == 1], 
                              xy[cp == 0], 
                              alternative = alt, 
                              exact = FALSE)  ## test stat positive if "1 > 0"
            res[rr, cc] <- tt$p.value
        }
    }
    expect_equal(fwt$p.value, res)
})


test_that("correctness of rowWilcoxonTests alternative", {
    p <- 30
    n <- 1000
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol=n)
    cls <- rbinom(n, 1, 1/2)
    
    ## adding some signal
    i1 <- which(cls==1)
    idx_greater <- 1:10
    idx_less <- 11:20
    idx_two.sided <- c(idx_greater, idx_less)
    idx_null <- setdiff(1:p, idx_two.sided)
    mat[idx_greater, i1] <- 1 + mat[idx_greater, i1]  
    mat[idx_less, i1] <- -1 + mat[idx_less, i1]  
    
    alt <- "two.sided"    
    pval <- rowWilcoxonTestsV1(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mt, m0)
    expect_lt(ml, m0)
    expect_lt(mg, m0)
    
    alt <- "greater"    
    pval <- rowWilcoxonTestsV1(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mg, m0)
    expect_gt(ml, m0)
    
    alt <- "less"    
    pval <- rowWilcoxonTestsV1(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_gt(mg, m0)
    expect_lt(ml, m0)
})


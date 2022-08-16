context("Calculation of Welch test statistics and p-values")

p <- 53
n <- 35
## null distribution
mat <- matrix(rnorm(p*n), ncol=n)
cls <- rep(c(0, 1), times=c(10, n-10))

test_that("rowWelchTests <=> welch.test", {
    alt <- sample(c("two.sided", "greater", "less"), 1)
    
    # Ordinary Welch t-test (reasonably fast for small number of tests)
    wt <- apply(mat, 1, FUN = function(x) {
        tt <- t.test(x[cls == 1], x[cls == 0], alternative = alt)  ## test stat positive if "1 > 0"
        c(tt$p.value, tt$statistic, tt$parameter)
    })
    dwt <- as.data.frame(t(wt))
    names(dwt) <- c("p.value", "statistic", "parameter")
    
    fwt1 <- rowWelchTests1(mat, categ = cls, alternative = alt)
    fwt <- rowWelchTests(mat, categ = cls, alternative = alt)
    expect_equal(fwt1, fwt)

    expect_equal(fwt$p.value, dwt$p.value)
    expect_equal(fwt$statistic, dwt$statistic)
    expect_equal(fwt$parameter, dwt$parameter)
})

test_that("rowWelchTests <=> welch.test (several contrasts at a time)", {
    alt <- sample(c("two.sided", "greater", "less"), 1)
    mat_small <- head(mat, 20)
    
    cls_perm <- replicate(10, sample(cls))
    fwt <- rowWelchTests(mat_small, categ = cls_perm, alternative = alt)
    
    # Ordinary Welch t-test (reasonably fast for small number of tests)
    nr <- nrow(mat_small)
    nc <- ncol(cls_perm)
    res <- matrix(NA_real_, nr, nc)
    for (rr in seq_len(nr)) {
        xy <- mat_small[rr, ]
        for (cc in seq_len(nc)) {
            cp <- cls_perm[, cc]
            tt <- t.test(xy[cp == 1], 
                         xy[cp == 0], 
                         alternative = alt)  ## test stat positive if "1 > 0"
            res[rr, cc] <- tt$p.value
        }
    }
        
    expect_equivalent(fwt$p.value, res)
})

test_that("correctness of rowWelchTests alternative", {
    p <- 100
    n <- 54
    ## null distribution
    mat <- matrix(rnorm(p*n), ncol=n)
    cls <- rep(c(0, 1), times = c(n/2, n-n/2))
    
    ## adding some signal
    i1 <- which(cls == 1)
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

    replicate(10, {
        x <- rnorm(1234)
        y <- rnorm(2345)
        target <- t.test(x, y)
        swt <- suffWelchTests(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
        swt1 <- suffWelchTests1(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
        expect_identical(swt, swt1)

        expect_equivalent(swt$statistic, target$statistic)
        expect_equal(swt$p.value, target$p.value)
        expect_equivalent(swt$parameter, target$parameter)
        expect_equivalent(swt$estimate, -diff(target$estimate))
    })
})


test_that("Sanity checks of 'categCheck' throw errors when expected to", {
    n = 10
    
    vec = rep(1, 5)
    expect_error(categCheck(vec, n), 
                 "vec should be of length 10, not5")
    
    categCheck(c(0,1,0,1), 4)
    
    # expect_error(categCheck(c(1,2,1,2), 4),
    #              "'c1212' should consist only of '0' and '1'! Or disctinct values (for a covariate).")

    categCheck(c(1,2,3,4), 4)
})

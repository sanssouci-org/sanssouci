context("Calculation of Binomial test statistics and p-values")

test_that("rowBinomialTests <=> binom.test", {
    p <- 250
    n0 <- 40; n1 <- 60
    mat0 <- matrix(rbinom(p*n0, size = 1, prob = 0.05), ncol = n0)
    mat1 <- matrix(rbinom(p*n1, size = 1, prob = 0.02), ncol = n1)
    mat <- cbind(mat0, mat1)
    cls <- rep(c(0, 1), times = c(n0, n1))

    for (alt in c("two.sided", "greater", "less")) {
        if (alt == "two.sided") {
            expect_warning(fbt <- rowBinomialTests(mat = mat, categ = cls, alternative = alt))
        } else {
            fbt <- rowBinomialTests(mat = mat, categ = cls, alternative = alt)
        }
        
        # Ordinary Binomial test (reasonably fast for 200 tests)
        pbt <- t(sapply(1:p, FUN=function(ii) {
          x1 <- mat[ii, cls == 1]
          x0 <- mat[ii, cls == 0]
          bt <- binom.test(sum(x1), length(x1), mean(x0), alternative = alt)
          c(statistic = bt[["statistic"]], p.value = bt[["p.value"]])
        }))

        pbt <- as.data.frame(pbt)
        names(pbt) <- c("statistic", "p.value")
        expect_equal(fbt$p.value, pbt$p.value, tolerance = 1e-12)
        expect_equal(fbt$statistic, pbt$statistic, tolerance = 1e-12)
    }
})


test_that("correctness of rowBinomialTests alternative", {
    
    p <- 123
    n0 <- 600; n1 <- 400
    p0 <- 0.05
    mat0 <- matrix(rbinom(p*n0, size = 1, prob = p0), ncol = n0)

    idx_null <- 1:round(p/3)
    idx_greater <- (round(p/3)+1):round(2*p/3)
    idx_less <- (round(2*p/3)+1):p
    idx_two.sided <- c(idx_greater, idx_less)
    
    pn <- length(idx_null)
    mat1.0 <- matrix(rbinom(pn*n1, size = 1, prob = p0), ncol = n1)

    pg <- length(idx_greater)
    mat1.g <- matrix(rbinom(pg*n1, size = 1, prob = 2*p0), ncol = n1)

    pl <- length(idx_less)
    mat1.l <- matrix(rbinom(pl*n1, size = 1, prob = p0/2), ncol = n1)

    mat1 <- rbind(mat1.0, mat1.g, mat1.l)
    mat <- cbind(mat0, mat1)
    cls <- rep(c(0, 1), times = c(n0, n1))
    
    alt <- "two.sided"
    expect_warning(pval <- rowBinomialTests(mat, categ = cls, alternative = alt)$p.value)
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mt, m0)
    expect_lt(ml, m0)
    expect_lt(mg, m0)

    alt <- "greater"    
    pval <- rowBinomialTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_lt(mg, m0)
    expect_gt(ml, m0)
    
    alt <- "less"    
    pval <- rowBinomialTests(mat, categ = cls, alternative = alt)$p.value
    m0 <- mean(pval[idx_null])
    mg <- mean(pval[idx_greater])
    ml <- mean(pval[idx_less])
    mt <- mean(pval[idx_two.sided])
    expect_gt(mg, m0)
    expect_lt(ml, m0)
})


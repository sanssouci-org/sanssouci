context("Upper bound of the number of false positives in a selection")

m <- 123
sim <- gaussianSamples(m = m, rho = 0.2, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
X <- sim$X
categ <- sim$categ
thr <- SimesThresholdFamily(m)(alpha = 0.2)
tests <- rowWelchTests(X, categ)
pval <- sort(tests$p.value)

test_that("Upper bound on the number of false positives", {
    expect_error(maxFP(pval, rev(thr)))  ## thr is not sorted descendingly
    
    best <- head(pval)
    worse <- tail(pval)
    b1 <- maxFP(best, thr)
    b2 <- maxFP(c(best, worse), thr)
    b3 <- maxFP(worse, thr)
    expect_lt(b1, b2)
    expect_lt(b1, b3)
    expect_lte(b3, b2)
    
    ub <- curveMaxFP(p.values = pval, thr)
    expect_equal(ub[length(best)], b1)
    
    M0 <- maxFP(pval, thr)
    sM0 <- maxFP(pval, thr)
    expect_equal(M0, sM0)
    expect_equal(M0, ub[m])
})

test_that("curveMaxFP", {
    expect_error(curveMaxFP(pval, rev(thr)))  
    expect_error(curveMaxFP(rev(pval), rev(thr)))
    
    ub <- curveMaxFP(pval, thr)
    expect_length(ub, m)
    expect_identical(ub, sort(ub))
    for (kk in seq_len(m)) {
        bk <- maxFP(pval[seq_len(kk)], thr)
        expect_equal(bk, ub[kk])
    }
})

test_that("flavors of curveMaxFP", {
    ub <- curveMaxFP(pval, thr)
    ub16 <- curveMaxFP(pval, thr, flavor = "BNR2014")
    ub14 <- curveMaxFP(pval, thr, flavor = "BNR2014") 
    ub06 <- curveMaxFP(pval, thr, flavor = "Mein2006") 
    expect_identical(ub, ub16)
    expect_identical(ub, ub14)
    expect_identical(ub, ub06)  ## not sure this is always true?
})

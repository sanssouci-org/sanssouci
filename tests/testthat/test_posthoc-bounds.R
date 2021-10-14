context("Upper bound of the number of false positives in a selection")

m <- 123
sim <- gaussianSamples(m = m, rho = 0.2, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
X <- sim$X
categ <- sim$categ
thr <- t_linear(0.2, seq_len(m), m)
tests <- rowWelchTests(X, categ)
pval <- sort(tests$p.value)

test_that("maxFP: Upper bound on the number of false positives", {
    expect_error(maxFP_low(pval, rev(thr)))  ## thr is not sorted descendingly
    
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
    
    expect_equal(maxFP(numeric(0L), thr), 0)
    expect_equal(curveMaxFP(numeric(0L), thr), numeric(0))
    
    expect_equal(curveMaxFP(numeric(0L), thr), numeric(0))

    bounds <- posthoc_bound(p.values = pval, S = 1:m, thr = thr, lab = "method")
    bounds <- posthoc_bound(p.values = pval, S = 1, thr = thr, lab = "method")
    bounds <- posthoc_bound(p.values = pval, S = numeric(0), thr = thr, lab = "method")
    expect_equal(nrow(bounds), 0)
})


test_that("curveMaxFP", {
    expect_error(curveMaxFP_old(pval, rev(thr)))  
    expect_error(curveMaxFP_old(rev(pval), rev(thr)))
    
    ub <- curveMaxFP(pval, thr)
    expect_length(ub, m)
    expect_identical(ub, sort(ub))
    for (kk in seq_len(m)) {
        bk <- maxFP(pval[seq_len(kk)], thr)
        expect_equal(bk, ub[kk])
    }
})

test_that("flavors of curveMaxFP", {
    ub <- curveMaxFP_ECN(pval, thr)
    ub16 <- curveMaxFP_old(pval, thr, flavor = "BNR2014")
    ub14 <- curveMaxFP_old(pval, thr, flavor = "BNR2014") 
    ub06 <- curveMaxFP_old(pval, thr, flavor = "Mein2006") 
    expect_identical(ub, ub16)
    expect_identical(ub, ub14)
    expect_identical(ub, ub06)  ## not sure this is always true?
})

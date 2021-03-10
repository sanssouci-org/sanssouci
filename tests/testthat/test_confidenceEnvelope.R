context("Confidence envelopes on the false positives")


test_that("Consistency of output of 'confCurveFromFam'", {
    m <- 123
    alpha <- 0.1
    sim <- gaussianSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    dex <- rowWelchTests(sim$X, sim$categ)
    conf_bound <- confCurveFromFam(dex$p.value, refFamily = "Simes", param = alpha, what = c("FP", "TP", "FDP", "TDP"))
    
    ## FP and TP
    FP <- subset(conf_bound, stat == "FP")$bound
    x <- subset(conf_bound, stat == "FP")$x
    TP <- subset(conf_bound, stat == "TP")$bound
    xx <- subset(conf_bound, stat == "TP")$x
    expect_equal(x, xx)
    expect_equal(FP + TP, x)
    
    ## FP in [0, m]?
    expect_gte(min(FP), 0)
    expect_lte(max(FP), m)
    
    ## FDP and TDP
    FDP <- subset(conf_bound, stat == "FDP")$bound
    x <- subset(conf_bound, stat == "FP")$x
    TDP <- subset(conf_bound, stat == "TDP")$bound
    xx <- subset(conf_bound, stat == "TDP")$x
    expect_equal(x, xx)
    expect_equal(FDP + TDP, rep(1, m))
    
    ## FDP in [0, 1]?
    expect_gte(min(FDP), 0)
    expect_lte(max(FDP), 1)
})

test_that("Consistency between outputs of 'fit', 'bound' and 'confCurveFromFam'", {
    m <- 123
    alpha <- 0.1
    obj <- SansSouciSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    cal <- fit(obj, alpha = alpha, B = 1e2, refFamily = "Simes")    
    cb0 <- bound(cal, all=TRUE)
    cb1 <- cal$output$conf_bound
    expect_identical(cb0, cb1)
    
    # using 'confCurveFromFam'
    lambda <- cal$output$lambda
    cb2 <- confCurveFromFam(pValues(cal), "Simes", param = lambda)
    expect_identical(cb0[, c("stat", "bound")], cb2[, c("stat", "bound")])
    
    # using 'confCurveFromFam' and 'fit(..., B=0)' to bypass calibration
    cb0 <- confCurveFromFam(pValues(cal), "Simes", param = alpha)
    cal0 <- fit(obj, alpha = alpha, B = 0, refFamily = "Simes")    
    cb00 <- bound(cal0, all = TRUE)
    expect_identical(cb0[, c("stat", "bound")], cb00[, c("stat", "bound")])
})

test_that("Ordering of confidence curves with/without calibration", {
    m <- 1230
    alpha <- 0.1
    obj <- SansSouciSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    
    cal <- fit(obj, alpha = alpha, B = 1e2, refFamily = "Simes")    
    cb1 <- bound(cal, all=TRUE)
    
    TP1 <- subset(cb1, stat == "TP")$bound
    TP0 <- subset(cb0, stat == "TP" )$bound
    diff <- TP1 - TP0
    expect_lte(min(diff), 0)
})

test_that("Vanilla test for 'plotConfCurve'", {
    sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    rwt <- rowWelchTests(sim$X, categ = sim$categ, alternative = "greater")
    
    ce <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = 0.1)
    p <- plotConfCurve(ce, xmax = 200)
    expect_true(inherits(p, "ggplot"))
})
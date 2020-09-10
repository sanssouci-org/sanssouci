context("Confidence envelopes on the false positives")


test_that("Consistency of output of 'confidenceEnvelope'", {
    m <- 123
    alpha <- 0.1
    sim <- gaussianSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    dex <- rowWelchTests(sim$X, sim$categ)
    conf_env <- confidenceEnvelope(dex$p.value, refFamily = "Simes", param = alpha, what = c("FP", "TP", "FDP", "TDP"))
    
    ## FP and TP
    FP <- subset(conf_env, stat == "FP")$bound
    x <- subset(conf_env, stat == "FP")$x
    TP <- subset(conf_env, stat == "TP")$bound
    xx <- subset(conf_env, stat == "TP")$x
    expect_equal(x, xx)
    expect_equal(FP + TP, x)
    
    ## FP in [0, m]?
    expect_gte(min(FP), 0)
    expect_lte(max(FP), m)
    
    ## FDP and TDP
    FDP <- subset(conf_env, stat == "FDP")$bound
    x <- subset(conf_env, stat == "FP")$x
    TDP <- subset(conf_env, stat == "TDP")$bound
    xx <- subset(conf_env, stat == "TDP")$x
    expect_equal(x, xx)
    expect_equal(FDP + TDP, rep(1, m))
    
    ## FDP in [0, 1]?
    expect_gte(min(FDP), 0)
    expect_lte(max(FDP), 1)
})

test_that("Ordering of confidence envelopes with/without calibration", {
    m <- 123
    alpha <- 0.1
    sim <- gaussianSamples(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    cal <- calibrateJER(sim$X, sim$categ, B = 1e3, alpha = alpha, refFamily = "Simes")    
    conf_env <- confidenceEnvelope(cal$p.value, "Simes", param = alpha)
    TP_Simes <- subset(cal$conf_env, stat == "TP")$bound
    TP <- subset(conf_env, stat == "TP" )$bound
    diff <- TP_Simes - TP
    expect_lte(min(diff), 0)
})

test_that("Vanilla test for 'plotConfidenceEnvelope'", {
    sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
    rwt <- rowWelchTests(sim$X, categ = sim$categ, alternative = "greater")
    
    ce <- confidenceEnvelope(rwt$p.value, refFamily = "Simes", param = 0.1)
    p <- plotConfidenceEnvelope(ce, xmax = 200)
    expect_true(inherits(p, "ggplot"))
})
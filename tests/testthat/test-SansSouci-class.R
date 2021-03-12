test_that("Correctness of the constructor of SansSouci", {
    m <- 54
    n <- 132
    set.seed(0xBEEF)
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                             pi0 = 0.8, SNR = 0, prob = 0.4)
    expect_s3_class(obj, "SansSouci")
    set.seed(0xBEEF)
    sim <- gaussianSamples(m = m, rho = 0, n = n, 
                           pi0 = 0.8, SNR = 0, prob = 0.4)
    obj1 <- SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H)
    expect_s3_class(obj1, "SansSouci")
    expect_identical(obj1, obj)
    
    expect_equal(names(obj), c("input", "parameters", "output"))
    expect_null(obj$parameters)
    expect_null(obj$output)
    
    Y <- obj$input$Y
    expect_equal(nHyp(obj), m)
    expect_equal(obj$input$m, m)
    expect_equal(nrow(Y), m)
    
    expect_equal(nObs(obj), n)
    expect_equal(ncol(Y), n)
    
    expect_error(SansSouci(Y = sim$X, groups = sim$categ+1),
                 "'groups' should consist only of '0' and '1'!")
    expect_error(SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H+1),
                 "'truth' should consist only of '0' and '1'!")
})


test_that("Correctness of elements of fitted  'SansSouci' object", {
    m <- 54
    n <- 132
    
    alpha <- 0.1
    B <- 25
    K <- m/2
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                        pi0 = 0.8, SNR = 0, prob = 0.4)
    alt <- "greater"
    fam <- "Beta"
    res <- fit(obj, alpha = alpha, B = B, K = K, alternative = alt, family = fam)
    expect_s3_class(res, "SansSouci")
    expect_identical(names(res), names(obj))
    params <- res$parameters
    names(params)
    expect_equal(params$alpha, alpha)
    expect_equal(params$B, B)
    expect_equal(params$alternative, alt)
    expect_equal(params$family, fam)
    expect_equal(params$K, K)
    
})


test_that("'fit.SansSouci' reproduces the results of 'calibrateJER'", {
    m <- 54
    n <- 132
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                         pi0 = 0.8, SNR = 0, prob = 0.4)
    Y <- obj$input$Y
    groups <- obj$input$groups
    alpha <- 0.07
    B <- 123
    K <- m/2
    alt <- "greater"
    fam <- "Beta"
    configs <- expand.grid(alternative = c("two.sided", "less", 
                                           "greater"), 
                           family = c("Simes", "Beta"), 
                           stringsAsFactors = FALSE)
    cc <- sample(nrow(configs), 1)   ## perform just one at random at each execution
    print(cc)
    alt <- configs[cc, "alternative"]
    fam <- configs[cc, "family"]
    set.seed(20210311)
    res <- fit(obj, alpha = alpha, B = B, K = K, 
               alternative = alt, family = fam)
    set.seed(20210311)
    cal <- calibrateJER(Y, groups, B = B, alpha = alpha, K = K, 
                        alternative = alt, refFamily = fam)
    expect_identical(cal, res$output)
})

test_that("Consistency of output of 'bound.SansSouci'", {
    m <- 54
    n <- 132
    
    alpha <- 0.1
    B <- 25
    K <- m/2
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                        pi0 = 0.8, SNR = 3, prob = 0.4)
    alt <- "greater"
    fam <- "Beta"
    res <- fit(obj, alpha = alpha, B = B, K = K, alternative = alt, family = fam)
    res <- fit(obj, alpha = alpha, B = B)
    what0 <- c("FP", "TP", "FDP", "TDP")
    
    # 'all=FALSE' => return a vector
    b <- predict(res, what = what0)
    expect_type(b, "double")
    expect_identical(names(b), what0)
    expect_equal(b[["FDP"]] + b[["TDP"]], 1)
    expect_equal(b[["FP"]] + b[["TP"]], m)
    
    # 'all=FALSE' => return a data.frame
    b <- predict(res, what = what0, all = TRUE)
    expect_s3_class(b, "data.frame")
    expect_equal(nrow(b), m * length(what0))
    
    # strict subset
    S <- order(pValues(res))[seq_len(m-10)]
    bb <- predict(res, S, what = what0, all = TRUE)
    expect_s3_class(bb, "data.frame")
    expect_equal(nrow(bb), length(S) * length(what0))
    names(bb)
    
    expect_equivalent(subset(b, x <= length(S)), bb)
    
    ww <- which(bb$x == length(S))
    w <- which(b$x == length(S))
    expect_equivalent(b[w, ], bb[ww, ])
})    


test_that("'bound.SansSouci' reproduces the results of 'curveMaxFP'", {
    m <- 54
    n <- 132
    
    alpha <- 0.1
    B <- 25
    K <- m/2
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                        pi0 = 0.8, SNR = 3, prob = 0.4)
    alt <- "greater"
    fam <- "Beta"
    res <- fit(obj, alpha = alpha, B = B, K = K, alternative = alt, family = fam)
    what0 <- c("FP", "TP", "FDP", "TDP")
    
    pvals <- sort(pValues(res))
    FP <- sansSouci:::curveMaxFP(p.values = pvals, 
                                 thr = thresholds(res))
    FPb <- predict(res, what = "FP", all = TRUE)$bound
    expect_identical(FPb, FP)
})
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
    
    output <- res$output
    nms <- c("statistic", "parameter", "p.value", "estimate", "p0", "thr",
             "piv_stat", "lambda", "steps_down")
    expect_identical(names(output), nms)
    expect_length(output$statistic, m)
    expect_length(output$parameter, m)
    expect_length(output$p.value, m)
    expect_length(output$estimate, m)
    p0 <- output$p0
    expect_equal(nrow(p0), m)
    expect_equal(ncol(p0), B)
    expect_lte(max(p0), 1)
    expect_gte(min(p0), 0)
    expect_length(output$thr, K)
    expect_length(output$piv_stat, B)
    expect_gte(output$lambda, 0)
    expect_lte(output$lambda, 1)
    expect_gte(output$steps_down, 0)
    
    expect_error(fit(obj, alpha = "alpha"))
    expect_error(fit(obj, alpha = "alpha"))
    expect_error(fit(obj))
})


test_that("'fit.SansSouci' reproduces the results of 'calibrate'", {
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
    alt <- configs[cc, "alternative"]
    fam <- configs[cc, "family"]
    set.seed(20210311)
    res <- fit(obj, alpha = alpha, B = B, K = K, 
               alternative = alt, family = fam, max_steps_down = 0)
    set.seed(20210311)
    p0 <- get_perm(obj$input$Y, obj$input$groups, B, alternative = alt)$p.value
    expect_equal(p0, res$output$p0)
    t_inv <- ifelse(fam == "Simes", t_inv_linear,  t_inv_beta)
    t_ <- ifelse(fam == "Simes", t_linear,  t_beta)
    pivStat <- get_pivotal_stat(p0, m, t_inv, K = K)
    expect_equal(pivStat, res$output$piv_stat)
#    expect_equal(quantile(pivStat, alpha), res$output$lambda)
#    expect_identical(cal$thr,      reso$thr)
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


test_that("'predict.SansSouci' reproduces the results of 'curveMaxFP'", {
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
    FP <- sanssouci:::curveMaxFP(p.values = pvals, 
                                 thr = thresholds(res))
    FPb <- predict(res, what = "FP", all = TRUE)$bound
    expect_identical(FPb, FP)
})

test_that("Correctness of Oracle predictions", {
    m <- 54
    pi0 <- 0.5
    m0 <- m * pi0
    m1 <- m - m0
    n <- 132
    
    alpha <- 0.05
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                        pi0 = pi0, SNR = 10, prob = 0.5)
    res_oracle <- fit(obj, alpha = alpha, family = "Oracle")
    FP <- predict(res_oracle, what = "FP", all = TRUE)$bound
    expect_equal(FP, c(rep(0, m1), 1:m0))
    
    # random selection of m/2 hypotheses should contain at least 
    # one true and one false positive with overwhelming proba
    S <- sample(m, m/2)
    FP <- predict(res_oracle, S = S, what = "FP")
    expect_lt(FP, length(S))
    expect_gt(FP, 0)

    # selection of first 10 hyps in the order of p-values
    # should contain only signal 
    p_values <- pValues(res_oracle)
    S <- head(order(p_values), 10)
    FP <- predict(res_oracle, S = S, what = "FP")
    expect_equal(FP, 0)

    # selection of last 10 hyps in the order of p-values
    # should contain only noise
    p_values <- pValues(res_oracle)
    S <- tail(order(p_values), 10)
    FP <- predict(res_oracle, S = S, what = "FP")
    expect_equal(FP, length(S))
    
    FP <- predict(res_oracle, S = integer(0L), what = "FP")
    expect_equal(FP, 0)
})
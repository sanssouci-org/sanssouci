test_that("Correctness of the constructor of SansSouci", {
    m <- 54
    n <- 13
    set.seed(0xBEEF)
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                             pi0 = 0.8, SNR = 0, prob = 0.4)
    set.seed(0xBEEF)
    sim <- gaussianSamples(m = m, rho = 0, n = n, 
                           pi0 = 0.8, SNR = 0, prob = 0.4)
    obj1 <- SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H)
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
    n <- 13
    
    alpha <- 0.1
    B <- 25
    K <- m/2
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                        pi0 = 0.8, SNR = 0, prob = 0.4)
    alt <- "greater"
    fam <- "Beta"
    res <- fit(obj, alpha = alpha, B = B, K = K, alternative = alt, family = fam)
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
    n <- 13
    obj <- SansSouciSim(m = m, rho = 0, n = n, 
                         pi0 = 0.8, SNR = 0, prob = 0.4)
    Y <- obj$input$Y
    groups <- obj$input$groups
    
    alt <- "greater"
    fam <- "Beta"
    configs <- expand.grid(alternative = c("two.sided", "less", 
                                           "greater"), 
                           family = c("Simes", "Beta"), 
                           stringsAsFactors = FALSE)
    for (cc in seq_len(nrow(configs))) {
        alt <- configs[cc, "alternative"]
        fam <- configs[cc, "family"]
        set.seed(20210311)
        res <- fit(obj, alpha = alpha, B = B, K = K, 
                   alternative = alt, family = fam)
        set.seed(20210311)
        cal <- calibrateJER(Y, groups, B = B, alpha = alpha, K = K, 
                            alternative = alt, refFamily = fam)
        expect_identical(cal, res$output)
    }
})


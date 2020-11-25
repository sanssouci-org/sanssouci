context("JER calibration (modular version")

test_that("JER calibration and get_pivotal_stat* yield identical pivotal statistics", {
    m <- 130
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 10
    
    ## reference: 'calibrateJER' 
    ## (without step down and with dummy alpha)
    set.seed(0xBEEF)
    system.time(cal0 <- calibrateJER(X, categ, B, 
                                     maxStepsDown = 0, alpha = 1))
    set.seed(0xBEEF)
    system.time(cal0_B <- calibrateJER(X, categ, B, 
                                      refFamily = "Beta",
                                      maxStepsDown = 0, alpha = 1))
    
    set.seed(0xBEEF)
    
    tf_slow <- tempfile()
    Rprof(tf_slow)
    pivStat_slow <- get_pivotal_stat_slow(X, categ, B)
    set.seed(0xBEEF)
    pivStat_slow_B <- get_pivotal_stat_slow(X, categ, B, t_inv = t_inv_beta)
    Rprof(NULL)
    summaryRprof(tf_slow)
    
    expect_equal(cal0$pivStat, pivStat_slow)
    expect_equal(cal0_B$pivStat, pivStat_slow_B)
    
    set.seed(0xBEEF)
    tf <- tempfile()
    Rprof(tf)
    pivStat_fast <- get_pivotal_stat_fast(X, categ, B)
    set.seed(0xBEEF)
    pivStat_fast_B <- get_pivotal_stat_fast(X, categ, B, t_inv = t_inv_beta)
    Rprof(NULL)
    summaryRprof(tf)
    
    expect_equal(cal0$pivStat, pivStat_fast)
    expect_equal(cal0_B$pivStat, pivStat_fast_B)
    
    set.seed(0xBEEF)
    tf <- tempfile()
    Rprof(tf)
    p0 <- get_perm_p(X, categ, B)
    pivStat <- get_pivotal_stat(p0)
    pivStat_B <- get_pivotal_stat(p0, t_inv = t_inv_beta)  # here we don't need to reset the seed!
    Rprof(NULL)
    summaryRprof(tf)
    
    expect_equal(cal0$pivStat, pivStat)
    expect_equal(cal0_B$pivStat, pivStat_B)
})

test_that("get_pivotal_stat* and get_one_pivotal_stat* yield identical results", {
    m <- 130
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 10
    
    set.seed(0xBEEF)
    tf <- tempfile()
    Rprof(tf)
    system.time(pivStat <- get_pivotal_stat_fast(X, categ, B))
    Rprof(NULL)
    summaryRprof(tf)
    
    set.seed(0xBEEF)
    tf_slow <- tempfile(fileext = ".Rout")
    Rprof(tf_slow)
    pivStat1_slow <- replicate(B, 
                               get_one_pivotal_stat_slow(X, sample(categ)))
    Rprof(NULL)
    summaryRprof(tf_slow)
    
    set.seed(0xBEEF)
    tf_slow <- tempfile()
    Rprof(tf_slow)
    pivStat1 <- replicate(B,
                          get_one_pivotal_stat(X, sample(categ)))
    Rprof(NULL)
    summaryRprof(tf_slow)
    
    expect_equal(pivStat, pivStat1_slow)  
    
    m <- 130
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 10
    
    ## reference: 'calibrateJER'
    ## (without step down and with dummy alpha)
    set.seed(0xBEEF)
    cal0 <- calibrateJER(X, categ, B, maxStepsDown = 0, alpha = 1)
    
    ## proposal:
    set.seed(0xBEEF)
    pivStat <- replicate(B, get_one_pivotal_stat_slow(X, sample(categ)))
    
    max(abs(cal0$pivStat-pivStat))
})
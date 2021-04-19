context("Joint Error Rate calibration")

test_that("Consistency of 'get_perm'", {
    m <- 5
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)

    B <- 10
    set.seed(123)
    perm0 <- get_perm(X, categ, B, rowWelchTests)

    set.seed(123)
    null_groups <- replicate(B, sample(categ))
    perm <- rowWelchTests(X, null_groups)
    expect_identical(perm0$p.value, perm$p.value)
    expect_identical(perm0$statistic, perm$statistic)

    
    set.seed(123)
    perm0 <- get_perm(X, categ, B, rowWilcoxonTests)
    set.seed(123)
    perm1 <- get_perm(X, categ, B, rowWilcoxonTests1)
    expect_identical(perm0$p.value, perm1$p.value)
    expect_identical(perm0$statistic, perm1$statistic)

    perm <- rowWilcoxonTests(X, null_groups)
    expect_identical(perm0$p.value, perm$p.value)
    expect_identical(perm0$statistic, perm$statistic)
})

test_that("Vanilla tests for 'get_pivotal_stat'", {
    m <- 50
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    w1 <- which(categ==1)
    X [, w1] <- X [, w1] + 0
    B <- 100
    null_groups <- replicate(B, sample(categ))
    p0 <- rowWelchTests(X, null_groups)$p.value
    
    pivStat <- get_pivotal_stat(p0)
    expect_length(pivStat, B)
    expect_lte(max(pivStat), 1)
    expect_gte(min(pivStat), 0) 

    pivStat2 <- get_pivotal_stat(p0, K = m/10)
    expect_gte(min(pivStat2), 0)
    expect_gte(min(pivStat2 - pivStat), 0)
})

test_that("JER calibration and get_pivotal_stat yield identical pivotal statistics", {
    m <- 132
    n <- 54
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 11
    
    set.seed(0xBEEF)
    p0 <- get_perm(X, categ, B)$p.value
    pivStat <- get_pivotal_stat(p0, t_inv_linear)
    
    # reference: 'calibrateJER' 
    # (without step down and with dummy alpha)
    set.seed(0xBEEF)
    cal0 <- calibrateJER(X, categ, B, 
                         maxStepsDown = 0, alpha = 1)
    expect_equal(cal0$piv_stat, pivStat)
    
    # Beta template
    set.seed(0xBEEF)
    cal0 <- calibrateJER(X, categ, B, 
                         refFamily = "Beta",
                         maxStepsDown = 0, alpha = 1)
    pivStat <- get_pivotal_stat(p0, t_inv_beta)
    expect_equal(cal0$piv_stat, pivStat)
    
    # K < m
    K <- round(2*m/3)
    set.seed(0xBEEF)
    cal0 <- calibrateJER(X, categ, B, 
                         maxStepsDown = 0, alpha = 1, K = K)
    pivStat_m <- get_pivotal_stat(p0, t_inv_linear, m = m, K = K)
    pivStat <- get_pivotal_stat(p0, t_inv_linear, K = K)
    expect_equal(cal0$piv_stat, pivStat_m)
    expect_equal(cal0$piv_stat, pivStat)
})



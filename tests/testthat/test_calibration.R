context("Joint Error Rate calibration")

test_that("Consistency of 'get_randomized_p_values_two_sample'", {
    m <- 5
    n <- 45
    X <- matrix(rnorm(m*n, mean = 1), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)

    B <- 10
    set.seed(123)
    perm0 <- get_randomized_p_values_two_sample(X, categ, B, rowWelchTests)

    set.seed(123)
    null_groups <- replicate(B, sample(categ))
    perm <- rowWelchTests(X, null_groups)
    expect_identical(perm0, perm$p.value)

    
    set.seed(123)
    perm0 <- get_randomized_p_values_two_sample(X, categ, B, rowWilcoxonTests)
    set.seed(123)
    perm1 <- get_randomized_p_values_two_sample(X, categ, B, rowWilcoxonTests1V1)
    expect_identical(perm0, perm1)

    perm <- rowWilcoxonTests(X, null_groups)
    expect_identical(perm0, perm$p.value)
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
    
    pivStat <- get_pivotal_stat(p0, m)
    expect_length(pivStat, B)
    expect_lte(max(pivStat), 1)
    expect_gte(min(pivStat), 0) 

    pivStat2 <- get_pivotal_stat(p0, m, K = m/10)
    expect_gte(min(pivStat2), 0)
    expect_gte(min(pivStat2 - pivStat), 0)
})

test_that("JER calibration and get_pivotal_stat yield identical pivotal statistics", {
    m <- 132
    n <- 54
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 11
    
    ## TODO: tests all param combinations
    set.seed(0xBEEF)
    p0 <- get_randomized_p_values_two_sample(X, categ, B)
    pivStat <- get_pivotal_stat(p0, m, t_inv_linear)

    cal <- calibrate(p0, m, alpha = 1, family = "Linear")
    expect_equal(cal$piv_stat, pivStat)
})

test_that("Test of get_randomized_p_values_one_sample", {
  n <- 10
  m <- 50
  set.seed(42)
  X <- matrix(rnorm(m*n), ncol = n, nrow = m)
  B <- 100
  pvals <- get_randomized_p_values_one_sample(X, B=B)
  expect_equal(dim(pvals), c(m,B))
  expect_equal(pvals, apply(pvals, 2, sort))
  expect_gt(min(pvals), 1e-7) #0 ? 
  expect_lte(max(pvals), 1)
  expect_lt(sum(pvals<0.1), B*m*0.12) #why this one ? 
})



test_that("pivotal statistics are non-decreasing in the set of hypotheses used for calibration", {
    m <- 132
    n <- 54
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 111
    
    p0 <- get_randomized_p_values_two_sample(X, categ, B)
    piv_stat <- get_pivotal_stat(p0, m, t_inv_linear)
    
    sub <- sample(m, round(m/2))
    piv_stat_sub <- get_pivotal_stat(p0[sub, ], m, t_inv_linear)
    
    expect_gte(min(piv_stat_sub - piv_stat_sub), 0)
})

test_that("Test consistency between R and python implementations", {
  p_values <- t(read.csv("test_cal_small.csv", header=FALSE))
  cal <- calibrate(p0 = p_values, m = nrow(p_values), alpha = 0.1, family = "Simes")
  expect_equal(unname(cal$lambda), 0.1541749)
})



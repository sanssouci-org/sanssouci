context("JER calibration (modular version")

test_that("JER calibration and get_pivotal_stat yield identical pivotal statistics", {
    m <- 130
    n <- 45
    X <- matrix(rnorm(m*n), ncol = n, nrow = m)
    categ <- rbinom(n, 1, 0.4)
    B <- 10
    
    ## reference: 'calibrateJER' 
    ## (without step down and with dummy alpha)
    set.seed(0xBEEF)
    cal0 <- calibrateJER(X, categ, B, 
                                     maxStepsDown = 0, alpha = 1)
    set.seed(0xBEEF)
    cal0_B <- calibrateJER(X, categ, B, 
                           refFamily = "Beta",
                           maxStepsDown = 0, alpha = 1)

    set.seed(0xBEEF)
    p0 <- get_perm_p(X, categ, B)
    pivStat <- get_pivotal_stat(p0)
    pivStat_B <- get_pivotal_stat(p0, t_inv = t_inv_beta)  # here we don't need to reset the seed!

    expect_equal(cal0$pivStat, pivStat)
    expect_equal(cal0_B$pivStat, pivStat_B)
})

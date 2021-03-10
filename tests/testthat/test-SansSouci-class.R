test_that("Correctness of the constructor of SansSouci", {
    m <- 54
    n <- 13
    set.seed(0xBEEF)
    obj1 <- SansSouciSim(m = m, rho = 0, n = n, 
                             pi0 = 0.8, SNR = 0, prob = 0.4)
    set.seed(0xBEEF)
    sim <- gaussianSamples(m = m, rho = 0, n = n, 
                           pi0 = 0.8, SNR = 0, prob = 0.4)
    obj <- SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H)
    expect_identical(obj1, obj)
    expect_equal(nHyp(obj), m)
    expect_equal(nObs(obj), n)
    expect_error(SansSouci(Y = sim$X, groups = sim$categ+1),
                 "'groups' should consist only of '0' and '1'!")
    expect_error(SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H+1),
                 "'truth' should consist only of '0' and '1'!")
    
})

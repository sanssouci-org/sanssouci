context("Confidence envelopes on the false positives")

test_that("Basic checks for 'conf_curve_from_fam'", {
  m <- 52
  sim <- sim_gaussian_equi(
    m = m, rho = 0.5, n = 10, pi0 = 0.8, SNR = 3,
    prob = 0.5
  )
  rwt <- row_welch_tests(sim$X, sim$categ, alternative = "greater")
  p_vals <- rwt$p_value

  cb <- conf_curve_from_fam(p_vals, ref_family = "Simes", param = 0.1,
                            what = "TP")
  expect_equal(nrow(cb), m)
  expect_equal(cb$label, rep("Simes(0.1)", m))

  cb_lin <- conf_curve_from_fam(p_vals,
    ref_family = "Linear", param = 0.1,
    what = "TP"
  )
  expect_equal(cb$stat, cb_lin$stat)

  cb <- conf_curve_from_fam(p_vals,
    ref_family = "Simes", param = 0.1,
    what = c("TP", "FP")
  )
  expect_equal(nrow(cb), m * 2)

  cb <- conf_curve_from_fam(p_vals, ref_family = "Beta", param = 0.1,
                            what = "TP")
  expect_equal(nrow(cb), m)

  cb <- conf_curve_from_fam(p_vals,
    ref_family = "Oracle", param = rep(0, m),
    what = "TP"
  )
  expect_equal(nrow(cb), m)
  expect_error(conf_curve_from_fam(rwt$p_value,
    ref_family = "unknown_fam",
    param = 0.1
  ))
})

test_that("Consistency of output of 'conf_curve_from_fam'", {
  m <- 123
  alpha <- 0.1
  sim <- sim_gaussian_equi(
    m = m, rho = 0.5, n = 100, pi0 = 0.8,
    SNR = 3, prob = 0.5
  )
  dex <- row_welch_tests(sim$X, sim$categ)
  conf_bound <- conf_curve_from_fam(dex$p_value,
    ref_family = "Simes",
    param = alpha,
    what = c("FP", "TP", "FDP", "TDP")
  )

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

test_that(
  "Consistency between outputs of 'fit', 'bound' and 'conf_curve_from_fam'",
  {
    m <- 123
    alpha <- 0.1
    obj <- SansSouciSim(
      m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
      prob = 0.5
    )
    cal <- fit(obj, alpha = alpha, B = 1e2, family = "Simes")
    cb0 <- predict(cal, all = TRUE)

    # using 'conf_curve_from_fam'
    lambda <- cal$output$lambda
    cb2 <- conf_curve_from_fam(get(p_values(cal), "Simes", param = lambda)
    expect_identical(cb0[, c("stat", "bound")], cb2[, c("stat", "bound")])

    # using 'conf_curve_from_fam' and 'fit(..., B=0)' to bypass calibration
    cb0 <- conf_curve_from_fam(get_p_values(cal), "Simes", param = alpha)
    cal0 <- fit(obj, alpha = alpha, B = 0, family = "Simes")
    cb00 <- predict(cal0, all = TRUE)
    expect_identical(cb0[, c("stat", "bound")], cb00[, c("stat", "bound")])
  }
)

test_that("Ordering of confidence curves with/without calibration", {
  m <- 1230
  alpha <- 0.1
  obj <- SansSouciSim(m = m, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)

  cal0 <- fit(obj, alpha = alpha, B = 0, family = "Simes")
  cb0 <- predict(cal0, all = TRUE)
  cal <- fit(obj, alpha = alpha, B = 1e2, family = "Simes")
  cb1 <- predict(cal, all = TRUE)

  TP1 <- subset(cb1, stat == "TP")$bound
  TP0 <- subset(cb0, stat == "TP")$bound
  diff <- TP1 - TP0
  expect_lte(min(diff), 0)
})

test_that("Vanilla test for 'plot_conf_curve'", {
  sim <- sim_gaussian_equi(
    m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
    prob = 0.5
  )
  rwt <- row_welch_tests(sim$X, categ = sim$categ, alternative = "greater")

  ce <- conf_curve_from_fam(rwt$p_value, ref_family = "Simes", param = 0.1)
  p <- plot_conf_curve(ce, xmax = 200)
  expect_true(inherits(p, "ggplot"))
})

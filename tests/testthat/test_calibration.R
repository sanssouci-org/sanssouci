context("Joint Error Rate calibration")

test_that("Consistency of 'perm_stat'", {
  m <- 5
  n <- 45
  X <- matrix(rnorm(m * n, mean = 1), ncol = n, nrow = m)
  categ <- rbinom(n, 1, 0.4)

  B <- 10
  set.seed(123)
  perm0 <- perm_stat(X, categ, B, row_welch_tests)

  set.seed(123)
  null_groups <- replicate(B, sample(categ))
  perm <- row_welch_tests(X, null_groups)
  expect_equal(perm0$p.value, perm$p.value)
  expect_equal(perm0$statistic, perm$statistic)

  set.seed(123)
  perm0 <- perm_stat(X, categ, B, row_wilcoxon_tests)
  set.seed(123)
  perm1 <- perm_stat(X, categ, B, row_wilcoxon_tests1_v1)
  expect_identical(perm0$p.value, perm1$p.value)
  expect_identical(perm0$statistic, perm1$statistic)

  perm <- row_wilcoxon_tests(X, null_groups)
  expect_identical(perm0$p.value, perm$p.value)
  expect_identical(perm0$statistic, perm$statistic)
})

test_that("Vanilla tests for 'pivotal_stat'", {
  m <- 50
  n <- 45
  X <- matrix(rnorm(m * n), ncol = n, nrow = m)
  categ <- rbinom(n, 1, 0.4)
  w1 <- which(categ == 1)
  X[, w1] <- X[, w1] + 0
  B <- 100
  null_groups <- replicate(B, sample(categ))
  p0 <- row_welch_tests(X, null_groups)$p.value

  piv_stat <- pivotal_stat(p0, m)
  expect_length(piv_stat, B)
  expect_lte(max(piv_stat), 1)
  expect_gte(min(piv_stat), 0)

  piv_stat2 <- pivotal_stat(p0, m, K = m / 10)
  expect_gte(min(piv_stat2), 0)
  expect_gte(min(piv_stat2 - piv_stat), 0)
})

test_that("JER calibration and pivotal_stat yield identical pivotal
          statistics", {
            m <- 132
            n <- 54
            X <- matrix(rnorm(m * n), ncol = n, nrow = m)
            categ <- rbinom(n, 1, 0.4)
            B <- 11

            ## TODO: tests all param combinations
            set.seed(0xBEEF)
            p0 <- perm_stat(X, categ, B)$p.value
            piv_stat <- pivotal_stat(p0, m, t_inv_linear)

            cal <- calibrate(p0, m, alpha = 1, family = "Linear")
            expect_equal(cal$piv_stat, piv_stat)
          })


test_that("pivotal statistics are non-decreasing in the set of hypotheses used
          for calibration", {
            m <- 132
            n <- 54
            X <- matrix(rnorm(m * n), ncol = n, nrow = m)
            categ <- rbinom(n, 1, 0.4)
            B <- 111

            p0 <- perm_stat(X, categ, B)$p.value
            piv_stat <- pivotal_stat(p0, m, t_inv_linear)

            sub <- sample(m, round(m / 2))
            piv_stat_sub <- pivotal_stat(p0[sub, ], m, t_inv_linear)

            expect_gte(min(piv_stat_sub - piv_stat_sub), 0)
          })

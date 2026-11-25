context("Calculation of Binomial test statistics and p-values")

test_that("row_binomial_tests <=> binom.test", {
  p <- 25
  n0 <- 40
  n1 <- 60
  mat0 <- matrix(rbinom(p * n0, size = 1, prob = 0.05), ncol = n0)
  mat1 <- matrix(rbinom(p * n1, size = 1, prob = 0.02), ncol = n1)
  mat <- cbind(mat0, mat1)
  cls <- rep(c(0, 1), times = c(n0, n1))

  alt <- sample(c("two.sided", "greater", "less"), 1)
  if (alt == "two.sided") {
    expect_warning(fbt <- row_binomial_tests(
      mat = mat, categ = cls,
      alternative = alt
    ))
  } else {
    fbt <- row_binomial_tests(mat = mat, categ = cls, alternative = alt)
  }

  # Ordinary Binomial test (reasonably fast for small number of tests)
  pbt <- t(sapply(1:p, FUN = function(ii) {
    x1 <- mat[ii, cls == 1]
    x0 <- mat[ii, cls == 0]
    bt <- binom.test(sum(x1), length(x1), mean(x0), alternative = alt)
    c(
      statistic = bt[["statistic"]],
      p.value = bt[["p.value"]],
      estimate = bt[["estimate"]]
    )
  }))
  pbt <- as.data.frame(pbt)
  names(pbt) <- c("statistic", "p.value", "estimate")
  expect_equal(fbt$p.value, pbt$p.value)
  expect_equal(fbt$statistic, pbt$statistic)
  py <- rowMeans(mat[, cls == 0])
  expect_equal(fbt$estimate, pbt$estimate - py)

  # matrix version
  cls_perm <- cbind(cls, replicate(3, sample(cls)))
  fbt_perm <- row_binomial_tests(
    mat = mat, categ = cls_perm,
    alternative = alt,
    warn = FALSE
  )

  expect_equal(fbt_perm$statistic[, 1], fbt$statistic)
  expect_equal(fbt_perm$p.value[, 1], fbt$p.value)
  expect_equal(fbt_perm$estimate[, 1], fbt$estimate)
  expect_equal(ncol(fbt_perm$p.value), ncol(cls_perm))
})


test_that("correctness of row_binomial_tests alternative", {
  p <- 123
  n0 <- 600
  n1 <- 400
  p0 <- 0.05
  mat0 <- matrix(rbinom(p * n0, size = 1, prob = p0), ncol = n0)

  idx_null <- 1:round(p / 3)
  idx_greater <- (round(p / 3) + 1):round(2 * p / 3)
  idx_less <- (round(2 * p / 3) + 1):p
  idx_two_sided <- c(idx_greater, idx_less)

  pn <- length(idx_null)
  mat1_null <- matrix(rbinom(pn * n1, size = 1, prob = p0), ncol = n1)

  pg <- length(idx_greater)
  mat1_greater <- matrix(rbinom(pg * n1, size = 1, prob = 2 * p0), ncol = n1)

  pl <- length(idx_less)
  mat1_less <- matrix(rbinom(pl * n1, size = 1, prob = p0 / 2), ncol = n1)

  mat1 <- rbind(mat1_null, mat1_greater, mat1_less)
  mat <- cbind(mat0, mat1)
  cls <- rep(c(0, 1), times = c(n0, n1))

  alt <- "two.sided"
  expect_warning({
    pval <- row_binomial_tests(mat, categ = cls, alternative = alt)$p.value
  })
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_lt(mt, m0)
  expect_lt(ml, m0)
  expect_lt(mg, m0)

  alt <- "greater"
  pval <- row_binomial_tests(mat, categ = cls, alternative = alt)$p.value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_lt(mg, m0)
  expect_gt(ml, m0)

  alt <- "less"
  pval <- row_binomial_tests(mat, categ = cls, alternative = alt)$p.value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_gt(mg, m0)
  expect_lt(ml, m0)
})

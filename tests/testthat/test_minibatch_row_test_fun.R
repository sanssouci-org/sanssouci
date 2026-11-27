context("Calculation of mini batch ")


test_that("Testing mini batch on small datasets", {
  alt <- sample(c("two.sided", "greater", "less"), 1)

  p <- 53
  n <- 35
  ## null distribution
  mat <- matrix(rnorm(p * n), ncol = n)
  groups <- rep(c(0, 1), times = c(10, n - 10))

  tests <- c(row_welch_tests, row_wilcoxon_tests)

  for (test in tests) {
    fwt <- test(mat, groups, alternative = alt)
    mb <- minibatch_row_test_fun(
      row_test_fun = test, Y = mat, groups = groups,
      alternative = alt, max_batch_size = 1e6
    )
    expect_equal(mb, as.matrix(fwt$p_value))
  }


  groups_perm <- replicate(10, sample(groups))

  for (test in tests) {
    fwt <- test(mat, groups_perm, alternative = alt)
    mb <- minibatch_row_test_fun(
      row_test_fun = test, Y = mat, groups = groups_perm,
      alternative = alt, max_batch_size = 1e3
    )
    expect_equal(mb, fwt$p_value)
  }

  for (test in tests) {
    fwt <- test(mat, groups_perm, alternative = alt)
    mb <- minibatch_row_test_fun(
      row_test_fun = test, Y = mat, groups = groups_perm,
      alternative = alt, max_batch_size = 1e2
    )
    expect_equal(mb, fwt$p_value)
  }
})

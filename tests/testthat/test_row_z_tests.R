test_that("consistency between row_z_tests and direct implementation", {
  p <- 15
  n <- 38
  mat <- matrix(rnorm(p * n, mean = 1), ncol = n)
  zt <- row_z_tests(mat, alternative = "greater")

  pvals <- apply(mat, 1, FUN = function(x) {
    stat <- sum(x) / sqrt(length(x))
    pnorm(stat, lower.tail = FALSE)
  })
  expect_equal(zt$p_values[, 1], pvals)

  # Sign flipping
  B <- 10
  eps <- replicate(B, rbinom(n, 1, 0.5) * 2 - 1) ## Rademacher
  zt_perm <- row_z_tests(mat, eps, alternative = "greater")

  stat_mat <- matrix(nrow = p, ncol = B)
  p_mat <- matrix(nrow = p, ncol = B)
  for (j in 1:p) {
    stat <- colSums(mat[j, ] * eps) / sqrt(n)
    stat_mat[j, ] <- stat
    p_mat[j, ] <- pnorm(stat, lower.tail = FALSE)
  }

  expect_equal(zt_perm$statistics, stat_mat)
  expect_equal(zt_perm$p_values, p_mat)
})

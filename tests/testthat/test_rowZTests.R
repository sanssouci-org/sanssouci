test_that("consistency between rowZTests and direct implementation", {
  p <- 15
  n <- 38
  mat <- matrix(rnorm(p * n, mean = 1), ncol = n)
  zt <- rowZTests(mat, alternative = "greater")

  pvals <- apply(mat, 1, FUN = function(x) {
     stat <- sum(x) / sqrt(length(x))
     pnorm(stat, lower.tail = FALSE)
  })
  expect_equal(zt$p.value[, 1], pvals)

  # Sign flipping
  B <- 10
  eps <- replicate(B, rbinom(n, 1, 0.5)*2 - 1)  ## Rademacher
  zt_perm <- rowZTests(mat, eps, alternative = "greater")

  stat_mat <- matrix(nrow = p, ncol = B)
  p_mat <- matrix(nrow = p, ncol = B)
  for (j in 1:p) {
    stat <- colSums(mat[j, ] * eps) / sqrt(n)
    stat_mat[j, ] <- stat
    p_mat[j, ] <- pnorm(stat, lower.tail = FALSE)
  } 
    
  expect_equal(zt_perm$statistic, stat_mat)
  expect_equal(zt_perm$p.value, p_mat)
})

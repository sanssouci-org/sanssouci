context("Calculation of randomization p-values")

m <- 513
n <- 38
mat <- matrix(rnorm(m * n), ncol = n)
categ <- rep(c(0, 1), times = c(27, n - 27))

test_that("Correctness of permutation p-values", {
  res <- test_by_randomization(X = mat, categ, B = 100, rand_p_value = TRUE)
  expect_equal(res$flavor, "perm")

  # consistency between parametric and randomization p-values
  expect_gt(cor(res$p_values, res$rand_p_value), 0.9)

  # Two different ways of getting permutation p-values
  null_stats <- cbind(res$null_statistics, res$statistics)
  sw <- sweep(abs(null_stats),
              MARGIN = 1,
              STATS = abs(res$statistics),
              FUN = ">=")
  p_values <- rowMeans(sw)
  expect_identical(res$rand_p_value, p_values)

  ## stochastic domination by Uniform distribution
  for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
    cs <- colSums(res$rand_null_p_values <= alpha)
    expect_true(mean(cs / m) <= alpha)
  }
})


test_that("Correctness of sign-flipping p-values", {
  colnames(mat) <- NULL
  res <- test_by_randomization(X = mat, B = 100, rand_p_value = TRUE)
  expect_equal(res$flavor, "flip")

  colnames(mat) <- seq_len(ncol(mat))
  res <- test_by_randomization(X = mat, B = 100, rand_p_value = TRUE)
  expect_equal(res$flavor, "flip")

  # Two different ways of getting sign flipping p-values give
  # identical results
  null_stats <- cbind(res$null_statistics, res$statistics)
  sw <- sweep(abs(null_stats),
              MARGIN = 1,
              STATS = abs(res$statistics),
              FUN = ">=")
  p_values <- rowMeans(sw)
  expect_identical(res$rand_p_values, p_values)

  # consistency between parametric and randomization p-values
  for (alt in c("two.sided", "less", "greater")) {
    res <- test_by_randomization(X = mat, B = 100,
                                 rand_p_value = TRUE,
                                 alternative = alt)
    expect_gt(cor(res$p_values, res$rand_p_values), 0.9)
  }

  ## stochastic domination by Uniform distribution
  for (alpha in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
    cs <- colSums(res$rand_null_p_values <= alpha)
    expect_true(mean(cs / m) <= alpha)
  }
})

test_that("Consistency of randomization p-values with different alternatives", {
  B <- 10
  for (alt in c("two.sided", "less", "greater")) {

    ## permutation
    categ <- rep(c(0, 1), times = c(27, n - 27))
    two_sided <- test_by_randomization(X = mat, categ, B = B,
                                       alternative = "two.sided")
    greater <- test_by_randomization(X = mat, categ, B = B,
                                     alternative = "greater")
    less <- test_by_randomization(X = mat, categ, B = B, alternative = "less")
    expect_equal(two_sided$p_values, 2 * pmin(greater$p_values, less$p_values))

    ## swapping the class labels
    categ_inv <- 1 - categ
    two_sided_inv <- test_by_randomization(X = mat, categ_inv, B = B,
                                           alternative = "two.sided")
    greater_inv <- test_by_randomization(X = mat, categ_inv, B = B,
                                         alternative = "greater")
    less_inv <- test_by_randomization(X = mat, categ_inv, B = B,
                                      alternative = "less")
    expect_equal(two_sided_inv$p_values, two_sided$p_values)
    expect_equal(greater_inv$p_values, less$p_values)
    expect_equal(less_inv$p_values, greater$p_values)
    
    ## randomization
    two_sided <- test_by_randomization(X = mat, categ_inv, B = B,
                                       alternative = "two.sided")
    greater <- test_by_randomization(X = mat, categ_inv, B = B,
                                     alternative = "greater")
    less <- test_by_randomization(X = mat, categ_inv, B = B,
                                  alternative = "less")
    expect_equal(two_sided$p_values, 2 * pmin(greater$p_values, less$p_values))
  }
})

test_that("Correctness of sanity checks", {
  m <- 5
  n <- 20
  mat <- matrix(rnorm(m * n), ncol = n)
  categ <- rep(c(0, 1), times = c(2, n - 2))
  expect_error(
    test_by_randomization(mat, categ, B = 1),
    ">2 observations from each sample required for two-sample tests"
  )
  categ[1] <- 12
  expect_error(test_by_randomization(mat, categ, B = 1),
               "Tests for more than 2 classes not implemented yet")

})

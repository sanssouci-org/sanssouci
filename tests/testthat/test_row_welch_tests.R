context("Calculation of Welch test statistics and p-values")

p <- 53
n <- 35
## null distribution
mat <- matrix(rnorm(p * n), ncol = n)
groups <- rep(c(0, 1), times = c(10, n - 10))

test_that("row_welch_tests <=> welch.test", {
  alt <- sample(c("two.sided", "greater", "less"), 1)

  # Ordinary Welch t-test (reasonably fast for small number of tests)
  wt <- apply(mat, 1, FUN = function(x) {
    tt <- t.test(x[groups == 1], x[groups == 0], alternative = alt)
    ## test stat positive if "1 > 0"
    c(tt$p.value, tt$statistic, tt$parameter)
  })
  dwt <- as.data.frame(t(wt))
  names(dwt) <- c("p.value", "statistic", "parameter")

  fwt1 <- row_welch_tests1(mat, groups = groups, alternative = alt)
  fwt <- row_welch_tests(mat, groups = groups, alternative = alt)
  expect_equal(fwt1, fwt)

  expect_equal(fwt$p_value, dwt$p.value)
  expect_equal(fwt$statistic, dwt$statistic)
  expect_equal(fwt$parameter, dwt$parameter)
})

test_that("row_welch_tests <=> welch.test (several contrasts at a time)", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  mat_small <- head(mat, 20)

  groups_perm <- replicate(10, sample(groups))
  fwt <- row_welch_tests(mat_small, groups = groups_perm, alternative = alt)

  # Ordinary Welch t-test (reasonably fast for small number of tests)
  nr <- nrow(mat_small)
  nc <- ncol(groups_perm)
  res <- matrix(NA_real_, nr, nc)
  for (rr in seq_len(nr)) {
    xy <- mat_small[rr, ]
    for (cc in seq_len(nc)) {
      cp <- groups_perm[, cc]
      tt <- t.test(xy[cp == 1],
        xy[cp == 0],
        alternative = alt
      ) ## test stat positive if "1 > 0"
      res[rr, cc] <- tt$p.value
    }
  }

  expect_equivalent(fwt$p_value, res)
})

test_that("correctness of row_welch_tests alternative", {
  p <- 100
  n <- 54
  ## null distribution
  mat <- matrix(rnorm(p * n), ncol = n)
  groups <- rep(c(0, 1), times = c(n / 2, n - n / 2))

  ## adding some signal
  i1 <- which(groups == 1)
  idx_greater <- 1:10
  idx_less <- 11:20
  idx_two_sided <- c(idx_greater, idx_less)
  idx_null <- setdiff(1:p, idx_two_sided)
  mat[idx_greater, i1] <- 1 + mat[idx_greater, i1]
  mat[idx_less, i1] <- -1 + mat[idx_less, i1]

  alt <- "two.sided"
  pval <- row_welch_tests(mat, groups = groups, alternative = alt)$p_value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_lt(mt, m0)
  expect_lt(ml, m0)
  expect_lt(mg, m0)

  alt <- "greater"
  pval <- row_welch_tests(mat, groups = groups, alternative = alt)$p_value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_lt(mg, m0)
  expect_gt(ml, m0)

  alt <- "less"
  pval <- row_welch_tests(mat, groups = groups, alternative = alt)$p_value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two_sided])
  expect_gt(mg, m0)
  expect_lt(ml, m0)
})


test_that("suffWelchTest gives identical results to welch.test", {
  replicate(10, {
    x <- rnorm(1234)
    y <- rnorm(2345)
    target <- t.test(x, y)
    swt <- suff_welch_tests(mean(x), mean(y), sd(x), sd(y),
                            length(x), length(y))
    swt1 <- suff_welch_tests1(
      mean(x), mean(y), sd(x), sd(y), length(x),
      length(y)
    )
    expect_identical(swt, swt1)

    expect_equivalent(swt$statistic, target$statistic)
    expect_equal(swt$p_value, target$p.value)
    expect_equivalent(swt$parameter, target$parameter)
    expect_equivalent(swt$estimate, -diff(target$estimate))
  })
})


test_that("Sanity checks of '.check_groups' throw errors when expected to", {
  n <- 10

  vec <- rep(1, 5)
  expect_error(
    .check_groups(vec, n),
    "vec should be of length 10, not5"
  )

  .check_groups(c(0, 1, 0, 1), 4)

  .check_groups(c(1, 2, 3, 4), 4)
})

test_that("Corner case with almost null within-group level variance", {
  file <- system.file("testdata", "issue-94_y.rds", package = "sanssouci")
  y <- readRDS(file)
  groups <- c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0)

  res_stat <- t.test(y ~ groups)
  res_sanssouci <- sanssouci::row_welch_tests(y, groups)
  names(res_sanssouci)
  expect_equivalent(res_stat$statistic, -res_sanssouci$statistic)
  # inconsistent direction??
  expect_equivalent(res_stat$p.value, res_sanssouci$p_value)
  expect_equivalent(res_stat$parameter, res_sanssouci$parameter)
})

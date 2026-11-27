context("Calculation of Pearson's correlation test statistics and p-values")

p <- 53
n <- 35
## null distribution
mat <- matrix(rnorm(p * n), ncol = n)
groups <- rnorm(n)

test_that("row_pearson_correlation_tests <=> cor.test", {
  alt <- sample(c("two.sided", "greater", "less"), 1)

  # Ordinary Welch t-test (reasonably fast for small number of tests)
  wt <- apply(mat, 1, FUN = function(x) {
    tt <- cor.test(x, groups, alternative = alt) ## test stat positive if "1 > 0"
    c(tt$p.value, tt$statistic)
  })
  dwt <- as.data.frame(t(wt))
  names(dwt) <- c("p.value", "statistic")

  fwt <- row_pearson_correlation_tests(mat, groups = groups, alternative = alt)

  expect_equal(fwt$p_value, dwt$p.value)
  expect_equal(fwt$statistic, dwt$statistic)
})


test_that(
  "row_pearson_correlation_tests <=> cor.test (several contrasts at a time)",
  {
    alt <- sample(c("two.sided", "greater", "less"), 1)
    mat_small <- mat
    B <- 10
    groups_perm <- replicate(B, sample(groups))
    fwt <- row_pearson_correlation_tests(mat_small,
      groups = groups_perm,
      alternative = alt
    )

    # Ordinary Pearson's correlation-test
    # (reasonably fast for small number of tests)
    nr <- nrow(mat_small)
    nc <- ncol(groups_perm)
    res <- matrix(NA_real_, nr, nc)
    stat <- matrix(NA_real_, nr, nc)
    est <- matrix(NA_real_, nr, nc)
    for (rr in seq_len(nr)) {
      xy <- mat_small[rr, ]
      for (cc in seq_len(nc)) {
        cp <- groups_perm[, cc]
        tt <- cor.test(xy,
          cp,
          alternative = alt
        ) ## test stat positive if "1 > 0"
        res[rr, cc] <- tt$p.value
        stat[rr, cc] <- tt$statistic
      }
    }

    expect_equivalent(fwt$p_value, res)
    expect_equivalent(fwt$statistic, stat)
  }
)

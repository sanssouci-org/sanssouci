context("Calculation of Welch test statistics and p-values")

p <- 53
n <- 35
## null distribution
mat <- matrix(rnorm(p*n), ncol=n)
cls <- rep(c(-1, 1), times=c(10, n-10))

test_that("rowTtestsOneSample <=> t.test", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  
  # Ordinary Welch t-test (reasonably fast for small number of tests)
  wt <- apply(mat, 1, FUN = function(x) {
    tt <- t.test(x*cls, alternative = alt)  ## test stat positive if "1 > 0"
    c(tt$p.value, tt$statistic, tt$parameter)
  })
  dwt <- as.data.frame(t(wt))
  names(dwt) <- c("p.value", "statistic", "parameter")
  
  fwt <- rowTtestsOneSample(mat, label = cls, alternative = alt)
  
  expect_equal(as.vector(fwt$p.value), dwt$p.value)
  expect_equal(as.vector(fwt$statistic), dwt$statistic)
})

test_that("rowTtestsOneSample <=> t.test (several contrasts at a time)", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  mat_small <- head(mat, 20)

  cls_perm <- replicate(10, sample(cls))
  fwt <- rowTtestsOneSample(mat_small, labels = cls_perm, alternative = alt)

  # Ordinary Welch t-test (reasonably fast for small number of tests)
  nr <- nrow(mat_small)
  nc <- ncol(cls_perm)
  res <- matrix(NA_real_, nr, nc)
  for (rr in seq_len(nr)) {
    xy <- mat_small[rr, ]
    for (cc in seq_len(nc)) {
      cp <- cls_perm[, cc]
      tt <- t.test(xy*cp,
                   alternative = alt)  ## test stat positive if "1 > 0"
      res[rr, cc] <- tt$p.value
    }
  }

  expect_equivalent(fwt$p.value, res)
})

test_that("correctness of rowTtestsOneSample alternative", {
  p <- 100
  n <- 54
  ## null distribution
  mat <- matrix(rnorm(p*n), ncol=n)
  cls <- rep(c(1, 1), times = c(n/2, n-n/2))

  ## adding some signal
  i1 <- which(cls == 1)
  idx_greater <- 1:10
  idx_less <- 11:20
  idx_two.sided <- c(idx_greater, idx_less)
  idx_null <- setdiff(1:p, idx_two.sided)
  mat[idx_greater, i1] <- 1 + mat[idx_greater, i1]
  mat[idx_less, i1] <- -1 + mat[idx_less, i1]

  alt <- "two.sided"
  pval <- rowTtestsOneSample(mat, labels = cls, alternative = alt)$p.value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two.sided])
  expect_lt(mt, m0)
  expect_lt(ml, m0)
  expect_lt(mg, m0)

  alt <- "greater"
  pval <- rowTtestsOneSample(mat, labels = cls, alternative = alt)$p.value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two.sided])
  expect_lt(mg, m0)
  expect_gt(ml, m0)

  alt <- "less"
  pval <- rowTtestsOneSample(mat, labels = cls, alternative = alt)$p.value
  m0 <- mean(pval[idx_null])
  mg <- mean(pval[idx_greater])
  ml <- mean(pval[idx_less])
  mt <- mean(pval[idx_two.sided])
  expect_gt(mg, m0)
  expect_lt(ml, m0)
})

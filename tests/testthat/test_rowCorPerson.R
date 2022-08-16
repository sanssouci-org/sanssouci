context("Calculation of Pearson's correlation test statistics and p-values")

p <- 53
n <- 35
## null distribution
mat <- matrix(rnorm(p*n), ncol=n)
cls <- rnorm(n)

test_that("rowCorPearson <=> cor.test", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  
  # Ordinary Welch t-test (reasonably fast for small number of tests)
  wt <- apply(mat, 1, FUN = function(x) {
    tt <- cor.test(x, cls, alternative = alt)  ## test stat positive if "1 > 0"
    c(tt$p.value, tt$statistic)
  })
  dwt <- as.data.frame(t(wt))
  names(dwt) <- c("p.value", "statistic")
  
  fwt <- rowCorPearson(mat, categ = cls, alternative = alt)
  
  expect_equal(fwt$p.value, dwt$p.value)
  expect_equal(fwt$statistic, dwt$statistic)
})


test_that("rowCorPearson <=> cor.test (several contrasts at a time)", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  mat_small <- mat
  B <- 10
  cls_perm <- replicate(B, sample(cls))
  fwt <- rowCorPearson(mat_small, categ = cls_perm, alternative = alt)
  
  # Ordinary Pearson's correlation-test (reasonably fast for small number of tests)
  nr <- nrow(mat_small)
  nc <- ncol(cls_perm)
  res <- matrix(NA_real_, nr, nc)
  stat <- matrix(NA_real_, nr, nc)
  est <- matrix(NA_real_, nr, nc)
  for (rr in seq_len(nr)) {
    xy <- mat_small[rr, ]
    for (cc in seq_len(nc)) {
      cp <- cls_perm[, cc]
      tt <- cor.test(xy, 
                   cp, 
                   alternative = alt)  ## test stat positive if "1 > 0"
      res[rr, cc] <- tt$p.value
      stat[rr,cc] <- tt$statistic
    }
  }
  
  expect_equivalent(fwt$p.value, res)
  expect_equivalent(fwt$statistic, stat)
  
})


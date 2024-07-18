context("Calculation of mini batch ")



test_that("Testing mini batch on small datasets", {
  alt <- sample(c("two.sided", "greater", "less"), 1)
  
  p <- 53
  n <- 35
  ## null distribution
  mat <- matrix(rnorm(p * n), ncol = n)
  cls <- rep(c(0, 1), times = c(10, n - 10))
  
  tests <- c(rowWelchTests, rowWilcoxonTests)
  
  for(test in tests){
    fwt <- test(mat, cls, alternative = alt)
    mb <- mini_batch_rowTestFUN(rowTestFUN = test, Y = mat, categ = cls,
                          alternative = alt, max_batch_size = 1e6)
    expect_equal(mb, as.matrix(fwt$p.value))
  }
  
  
  cls_perm <- replicate(10, sample(cls))
  
  for(test in tests){
    fwt <- test(mat, cls_perm, alternative = alt)
    mb <- mini_batch_rowTestFUN(rowTestFUN = test, Y = mat, categ = cls_perm,
                                alternative = alt, max_batch_size = 1e3)
    expect_equal(mb, fwt$p.value)
  }
  
  for(test in tests){
    fwt <- test(mat, cls_perm, alternative = alt)
    mb <- mini_batch_rowTestFUN(rowTestFUN = test, Y = mat, categ = cls_perm,
                                alternative = alt, max_batch_size = 1e2)
    expect_equal(mb, fwt$p.value)
  }
  
})

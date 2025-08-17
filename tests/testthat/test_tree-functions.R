context("Tree functions")

test_that("Vanilla test for tree functions", {
    m <- 6
    check_tree <- function(dd) {
        expect_true(class(dd) == "list")
        expect_length(dd, 2)
        expect_named(dd)
        expect_identical(names(dd), c("leaf_list", "C"))
    }

    meth <- sample(2, 1)
    dd <- dyadic.from.window.size(m, s = 2, method = meth)
    check_tree(dd)
    
    dd <- dyadic.from.height(m, H = 3, method = meth)
    check_tree(dd)
    
    dd <- dyadic.from.height(m, method = meth)
    check_tree(dd)
})

test_that("Vanilla test for 'zeta' functions", {
  m <- 100
  x <- rnorm(m, mean = c(rep(c(0, 2), each = 50)))
  pval <- 1 - pnorm(x)
  lambda <- 0.05
  
  expect_equal(zeta.trivial(pval, lambda), 
               length(x))
  
  expect_equal(zeta.HB(pval, lambda), 
               sum(p.adjust(pval, method = 'holm') >= lambda))
  
  expect_lte(zeta.DKWM(pval, lambda), m)
  expect_lte(zeta.HB(pval, lambda), zeta.kBonf(pval, lambda))
})

test_that("'curve.V.star.*' functions", {
  m <- 20
  C <- list(
    list(c(2, 5), c(8, 15), c(16, 19)),
    list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
    list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15), c(17, 17), c(18, 19)),
    list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
  )
  ZL <- list(
    c(4, 8, 4),
    c(3, 3, 4, 1, 3),
    c(2, 2, 1, 1, 2, 1, 2),
    c(1, 1, 1, 2, 1, 1)
  )
  leaf_list <- as.list(1:m)
  res_naive <- curve.V.star.forest.naive(1:m, C, ZL, leaf_list, pruning = FALSE)
  res_fast <- curve.V.star.forest.fast(1:m, C, ZL, leaf_list, pruning = FALSE)
  expect_equal(res_naive, res_fast)
  
  res_naive <- curve.V.star.forest.naive(1:m, C, ZL, leaf_list, pruning = TRUE)
  res_fast <- curve.V.star.forest.fast(1:m, C, ZL, leaf_list, pruning = TRUE)
  expect_equal(res_naive, res_fast)

  res_naive <- curve.V.star.forest.naive(1:m, C, ZL, leaf_list, pruning = TRUE, 
                                         delete.gaps = TRUE)
  res_fast <- curve.V.star.forest.fast(1:m, C, ZL, leaf_list, pruning = TRUE,
                                       delete.gaps = TRUE)
  expect_equal(res_naive, res_fast)
})



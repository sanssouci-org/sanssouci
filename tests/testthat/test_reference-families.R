test_that("Simes/Linear reference family", {
    m <- 10
    alpha <- 0.05
    k <- 5
    thr <- t_linear(alpha, 1:m, m)
    thr_inv <- t_inv_linear(thr, 1:m, m)
    expect_equal(thr_inv, rep(alpha, m))
    
})

test_that("Beta reference family", {
    m <- 10
    alpha <- 0.05
    k <- 5
    thr <- t_beta(alpha, 1:m, m)
    thr_inv <- t_inv_beta(thr, 1:m, m)
    expect_equal(thr_inv, rep(alpha, m))
})

test_that("check_ref_fam", {
  res <- check_ref_fam(Simes)
  expect_null(res)
  res <- check_ref_fam(Beta)
  expect_null(res)
})
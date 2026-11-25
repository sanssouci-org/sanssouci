context("Upper bound of the number of false positives in a selection")

m <- 123
sim <- sim_gaussian_equi(
  m = m, rho = 0.2, n = 100, pi0 = 0.8, SNR = 3,
  prob = 0.5
)
X <- sim$X
categ <- sim$categ
thr <- t_linear(0.2, seq_len(m), m)
tests <- row_welch_tests(X, categ)
pval <- sort(tests$p.value)
whats <- c("TP", "FP", "TDP", "FDP")

test_that("I/O of 'posthoc_bound'", {
  nb_stats <- ceiling(runif(1) * 4)
  what <- sample(whats, nb_stats)

  # subset of arbitrary size
  s <- sample(m, 1)
  S <- sample(m, s)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what
  )
  expect_equal(nrow(bounds), nb_stats)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what, all = TRUE
  )
  expect_equal(nrow(bounds), nb_stats * s)

  # subset of size 1
  s <- 1
  S <- sample(m, s)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what
  )
  expect_equal(nrow(bounds), nb_stats)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what, all = TRUE
  )
  expect_equal(nrow(bounds), nb_stats * s)

  # empty subset
  S <- integer(0L)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what
  )
  expect_equal(nrow(bounds), nb_stats)
  bounds <- posthoc_bound(
    p_values = pval, S = S, thr = thr, lab = "method",
    what = what, all = TRUE
  )
  expect_equal(nrow(bounds), 0)

  # not a subset
  S <- m + 1
  expect_error(posthoc_bound(p_values = pval, S = S, thr = thr, lab = "method"),
    regexp = "Argument 'S' should be a subset of indices between 1 and "
  )
  S <- 1 + 0:m
  expect_error(posthoc_bound(p_values = pval, S = S, thr = thr, lab = "method"),
    regexp = "Argument 'S' should be a subset of indices between 1 and "
  )

  # invalid 'p.values'
  pval_wrong <- pval
  pval_wrong[1] <- NA
  expect_error(posthoc_bound(p_values = pval_wrong, thr = thr, lab = "method"),
    regexp = "Argument 'p_values' should only contain elements between 0 and 1"
  )

  pval_wrong[1] <- -1
  expect_error(posthoc_bound(p_values = pval_wrong, thr = thr, lab = "method"),
    regexp = "Argument 'p_values' should only contain elements between 0 and 1"
  )

  pval_wrong[1] <- 2
  expect_error(posthoc_bound(p_values = pval_wrong, thr = thr, lab = "method"),
    regexp = "Argument 'p_values' should only contain elements between 0 and 1"
  )
})

test_that("max_fp and curve_max_fp are consistent", {
  expect_error(maxFP_low(pval, rev(thr))) ## thr is not sorted descendingly

  best <- head(pval)
  worse <- tail(pval)
  b1 <- max_fp(best, thr)
  b2 <- max_fp(c(best, worse), thr)
  b3 <- max_fp(worse, thr)
  expect_lt(b1, b2)
  expect_lt(b1, b3)
  expect_lte(b3, b2)

  ub <- curve_max_fp(p_values = pval, thr)
  expect_equal(ub[length(best)], b1)

  M0 <- max_fp(pval, thr)
  sM0 <- max_fp(pval, thr)
  expect_equal(M0, sM0)
  expect_equal(M0, ub[m])

  expect_equal(max_fp(numeric(0L), thr), 0)
  expect_equal(curve_max_fp(numeric(0L), thr), numeric(0))

  expect_equal(curve_max_fp(numeric(0L), thr), numeric(0))
})


test_that("curve_max_fp and max_fp are consistent for subsets", {
  ub <- curve_max_fp(pval, thr)
  expect_length(ub, m)
  expect_identical(ub, sort(ub))
  for (kk in seq_len(m)) {
    bk <- max_fp(pval[seq_len(kk)], thr)
    expect_equal(bk, ub[kk])
  }
})

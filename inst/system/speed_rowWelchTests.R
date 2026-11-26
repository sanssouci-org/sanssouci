p <- 1e5
n <- 383
mat <- matrix(rnorm(p * n), ncol = n)
cls <- rep(c(0, 1), times = c(274, n - 274))

system.time(fwt <- row_welch_tests(mat, categ = cls))
str(fwt)

# ordinary t.test:
system.time(pwt <- apply(mat, 1, FUN = function(x) {
  t.test(x[cls == 0], x[cls == 1])$p.value
}))

# compare:
stopifnot(sum(abs(fwt$p.value - pwt)) < 1e-10) ## same results

library(microbenchmark)
p <- 1e3
n <- 38
mat <- matrix(rnorm(p * n), ncol = n)
cls <- rep(c(0, 1), times = c(27, n - 27))

cls_perm <- replicate(100, sample(cls))
stats <- row_welch_tests(mat, categ = cls_perm)

# Compare to matrixTests (which can only handle one comparison, ie 'categ' must
# be a vector) from https://github.com/karoliskoncevicius/matrixTests
#
X <- matrix(rnorm(1e7), ncol = 1e1)
Y <- matrix(rnorm(1e7), ncol = 1e1)
system.time(res_matrixtests <- row_t_welch(X, Y)) # running time: 2.4 seconds

Z <- cbind(X, Y)
groups <- rep(c(1, 0), times = c(ncol(X), ncol(Y)))
system.time(res_sanssouci <- row_welch_tests(Z, groups))
# running time: <1 second

max(abs(res_matrixtests$statistic - res_sanssouci$statistic))
all.equal(res_matrixtests$statistic, res_sanssouci$statistic)

max(abs(res_matrixtests$pvalue - res_sanssouci$p.value))
all.equal(res_matrixtests$pvalue, res_sanssouci$p.value)

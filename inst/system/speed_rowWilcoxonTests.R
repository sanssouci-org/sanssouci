p <- 1e5
n <- 383
mat <- matrix(rnorm(p*n), ncol=n)
# mat <- matrix(rnbinom(n*p, size=1, prob = 0.1))
cls <- rep(c(0, 1), times=c(274, n-274))

system.time(fwt <- rowWilcoxonTests(mat, categ=cls))
str(fwt)

# ordinary t.test:
system.time(pwt <- apply(mat, 1, FUN=function(x) {
  wilcox.test(x[cls==0], x[cls==1])$p.value
}))

# compare:
stopifnot(sum(abs(fwt$p.value-pwt)) < 1e-10)  ## same results
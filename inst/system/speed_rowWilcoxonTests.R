p <- 1e5
n <- 383
# mat <- matrix(rnorm(p*n), ncol=n)
mat <- matrix(rnbinom(n*p, size=1, prob = 0.1), ncol=n)
cls <- (rep(c(0, 1), times=c(274, n-274)))

system.time(fwt <- rowWilcoxonTests(mat, categ=cls))
str(fwt)

# ordinary t.test:
system.time(pwt <- apply(mat, 1, FUN=function(x) {
  wilcox.test(x[cls==0], x[cls==1])$p.value
}))

# compare:
stopifnot(sum(abs(fwt$p.value-pwt)) < 1e-10)  ## same results

## comparison plot with microbenchmark
library(microbenchmark)

p <- 1e3
n <- 100
# mat <- matrix(rnorm(p*n), ncol=n)
mat <- matrix(rnbinom(n*p, size=1, prob = 0.1), ncol = n)
cls <- rep(c(0, 1), times=c(60, n-60))
mb <- microbenchmark(rowWilcoxonTest = {fwt <- rowWilcoxonTests(mat, categ=cls)},
                     wilcoxon.test = {pwt <- apply(mat, 1, FUN=function(x) {
                       wilcox.test(x[cls==0], x[cls==1],  exact = FALSE)$p.value
                     })},
                     times=10L)

ggplot2::autoplot(mb)

## comparison plot with microbenchmark
library(microbenchmark)

p <- 15
n <- 40
B = 10
# mat <- matrix(rnorm(p*n), ncol=n)
mat <- matrix(rnbinom(n*p, size=1, prob = 0.1), ncol = n)
cls <- rep(c(0, 1), times=c(20, n-20))
cls_mat <- replicate(B, sample(cls))
mb <- microbenchmark(rowWilcoxonTest = {fwt <- rowWilcoxonTests(mat, categ=cls_mat)},
                     wilcoxon.test = {
                       nr <- nrow(mat)
                       nc <- ncol(cls_mat)
                       res <- matrix(NA_real_, nr, nc)
                       for (rr in seq_len(nr)) {
                         xy <- mat[rr, ]
                         for (cc in seq_len(nc)) {
                           cp <- cls_mat[, cc]
                           tt <- wilcox.test(xy[cp == 1], 
                                             xy[cp == 0], 
                                             exact = FALSE)  ## test stat positive if "1 > 0"
                           res[rr, cc] <- tt$p.value
                         }
                       }
                     },
                     times=20L)

ggplot2::autoplot(mb)






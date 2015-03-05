library(sansSouci)

pi0 <- 1
m <- 1e2
p <- 0.5
n <- 100

nVio <- function(m, p, n, pi0, B=1e3, alpha=0.2, verbose=FALSE) {
    eps <- matrix(rnorm(m*n), m, n)
    y <- rbinom(n, 1, p)
    
    ## nulls and alternatives
    m0 <- round(m*pi0)
    m1 <- m-m0
    H0 <- 1:m0
    if (m0 < m) {
        H1 <- seq(from=m0+1, to=m, by=1)
    } else {
        H1 <- integer(0)
    }
  
    ## means under the total null: 0
    w <- wilcoxStat(eps, y, B=B)

    scoreMat <- w$stat0Mat
    o <- order(w$stat, decreasing=TRUE)
    stat <- w$stat[o]

    resJ <- getJointFWERThresholds(scoreMat, tau="kFWER", alpha=alpha, maxSteps=100)
    thr <- resJ$thr

    BB <- sapply(stat, function(x) sum(x<=thr))  ## Eqn (7) in Meinshausen (2006) (*not* Vbar)
    ##      R <- sapply(thr, function(x) sum(stat>=x))   ## Number of rejections
    R <- 1:m
    Sbar <- pmax(0, cummax(R-BB[R]))
    S <- cumsum(sapply(o, "%in%", H1))  ## number of TRUE discoveries among first rejections
    V <- cumsum(sapply(o, "%in%", H0))  ## number of FALSE discoveries among first rejections
    R <- S+V
    Vbar <- R-Sbar[R]                   ## (tighter) upper bound on number of FALSE discoveries 

    resS <- any(Sbar>S)
    resV <- any(Vbar<V)
    res <- any(stat>thr)
    if (verbose) {
        print(res)
        print(resV)
        print(resS)
    }
    stopifnot(resV==resS)
    if (pi0==1) {
        stopifnot(resS==res)
    }
    resS
}

run1 <- function(x, ...) nVio(m, p, n, pi0, ...)
run1()

nc <- 4
library(parallel)
system.time(res <- do.call(c, mclapply(1:10, run1, mc.cores=nc)))
mean(res>0)


run0 <- function(x, ...) nVio(m, p, n, pi0=1, ...)
run1 <- function(x, ...) nVio(m, p, n, pi0=0.99, ...)
run1(verbose=TRUE)

run10 <- function(x, ...) nVio(m, p, n, pi0=0.9, ...)
run10(verbose=TRUE)

if (FALSE) {
    print(system.time(res0 <- do.call(c, mclapply(1:1e3, run0, mc.cores=nc))))
    mean(res0>0)

    print(system.time(res1 <- do.call(c, mclapply(1:1e3, run1, mc.cores=nc))))
    mean(res1>0)

    system.time(res10 <- do.call(c, mclapply(1:1e2, run10, mc.cores=nc)))
    mean(res10>0)
}

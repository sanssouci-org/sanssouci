library(sansSouci)

pi0 <- 1
m <- 1e1
p <- 0.5
n <- 100

nVio <- function(m, p, n, pi0, B=1e3, alpha=0.2, stepDown=FALSE,
                 verbose=FALSE, browseVio=FALSE, plot=FALSE) {
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

    if (stepDown) {
        resJ <- stepDownControl(stat, scoreMat, tau="kFWER", alpha=alpha, verbose=verbose)
    } else {
        resJ <- getJointFWERThresholds(scoreMat, tau="kFWER", alpha=alpha)
    }
    thr <- resJ$thr

    BB <- sapply(stat, function(x) sum(x<=thr))  ## Eqn (7) in Meinshausen (2006) (*not* Vbar)
    ##      R <- sapply(thr, function(x) sum(stat>=x))   ## Number of rejections
    R <- 1:m
    Sbar <- pmax(0, cummax(R-BB[R]))
    S <- cumsum(sapply(o, "%in%", H1))  ## number of TRUE discoveries among first rejections
    V <- cumsum(sapply(o, "%in%", H0))  ## number of FALSE discoveries among first rejections
    R <- S+V
    Vbar <- R-Sbar[R]                   ## (tighter) upper bound on number of FALSE discoveries 

    if (plot) {
        par(lwd=2)
        plot(V, t='s')
        lines(Vbar, col=2, t='s')
        lines(1:m, col=3, t='s', lty=4)
    }
    
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
    if (browseVio) {
        browser()
    }
    c(res, resS)  ## 'res' assumes pi0=1, 'resS' does not
}

run1 <- function(x, ...) nVio(m, p, n, pi0, ...)
run1()

nc <- 4
library(parallel)
system.time(res <- do.call(cbind, mclapply(1:10, run1, mc.cores=nc)))
mean(res>0)


run0 <- function(x, ...) nVio(m, p, n, pi0=1, ...)
run1 <- function(x, ...) nVio(m, p, n, pi0=0.99, ...)
run1(verbose=TRUE)

run10 <- function(x, ...) nVio(m, p, n, pi0=0.9, ...)
run10(verbose=TRUE)

run20 <- function(x, ...) nVio(m, p, n, pi0=0.8, ...)
run20(plot=TRUE)


run20SD <- function(x, ...) nVio(m, p, n, pi0=0.8, stepDown=TRUE, verbose=TRUE, ...)
run20SD(plot=TRUE)

if (FALSE) {
    print(system.time(res0 <- do.call(cbind, mclapply(1:1e3, run0, mc.cores=nc))))
    rowMeans(res0)

    print(system.time(res1 <- do.call(cbind, mclapply(1:1e3, run1, mc.cores=nc))))
    rowMeans(res1)

    system.time(res10 <- do.call(cbind, mclapply(1:1e2, run10, mc.cores=nc)))
    rowMeans(res10)

    system.time(res20 <- do.call(cbind, mclapply(1:1e3, run20, mc.cores=nc)))
    rowMeans(res20)

    system.time(res20SD <- do.call(cbind, mclapply(1:1e3, run20SD, mc.cores=nc)))
    rowMeans(res20SD)

}


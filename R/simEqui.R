simEqui <- structure(function( ### Simulate Gaussian equi-correlated test statistics
        m,
        ### Number of hypotheses
        rho,
        ### Level of equi-correlation between pairs of variables
        B,
        ### Number of simulations
        pi0,
        ### Proportion of true null hypotheses
        SNR=1){
        ### Signal to noise ratio. Either a numeric value (a measure of distance between H0 and H1) or a vector of length \code{m*(1-pi0)}
    m0 <- round(m*pi0)
    m1 <- m-m0
    H <- rep(c(0, 1), times=c(m0, m1))
    H <- sample(H)
    H0 <- which(H==0)
    H1 <- which(H==1)
    if (length(SNR)>1) {
        stopifnot(length(SNR)==m1)
    }
    ## signals
    mu <- rep(0, m)
    mu[H1] <- SNR
    
    ## equi-correlated noise
    sim <- simulateGaussianEquiCorrelatedNulls(m, n=1+B, rho=rho)
    x <- mu + sim[, 1]
    Xb <- sim[, -1, drop=FALSE]
    list(
        x=x,
        ### A vector of length \eqn{m} test statistics
        X0=Xb,
        ### An \eqn{m x B} matrix of test statistics under the null
        H=H)
        ### A vector of length \eqn{m}, the status of each
        ### hypothesis:\describe{
        ### \item{0}{true null hypothesis}
        ### \item{1}{true alternative hypothesis}
        ### }
    }, ex=function() {
        m <- 123
        rho <- 0.2
        B <- 100
        pi0 <- 0.5
        
        sim <- simEqui(m, rho, B, pi0, SNR=1)
        X <- sim$X
        
        w <- wilcoxStat(X, y, B=B)
        scoreMat <- w$stat0Mat
        stat <- w$stat
        
        ## show test statistics
        pch <- 20
        colStat <- 1+sim$H
        plot(stat, col=colStat, main="Test statistics", pch=pch)
        legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
    })

##references<<Further results on controlling the false discovery
##proportion. Guo W., He L. Sarkar S. K. Theorem 3.8.

## G <- function(alpha, r, s, F) {
##   ar <- alpha[r]
##   ar1 <- alpha[r-1]
##   as <- alpha[s]
##   as1 <- alpha[s-1]
##   F(ar, as) - F(ar1, as) - F(ar, as1) + F(ar1, as1)
## }

mTilde <- function(i, n1, gamma) {
  pmin(i+n1, mStar2(i, gamma))
}

mStar <- function(i, gamma) {
  if (i==0) {
    res <- 0
  } else {
    idxs <- seq(from=1, to=n, by=1)
    res <- max(which(floor(gamma*idxs)+1 <= i))
  }
  res
}

mStar2 <- function(i, gamma) {
  res <- pmin(ceiling(i/gamma)-1, n)
  res[which(i==0)] <- 0
  res
}

mStar3 <- function(i, gamma) {
  res <- rep(n, length(i))
  ww <- which(i<n*gamma)
  if (length(ww)) {
    res[ww] <- ceiling(i[ww]/gamma)-1
  }
  res[which(i==0)] <- 0
  res
}

alphaPrime <- function(i, beta, gamma) {
  gammaI <- floor(gamma*i)+1
  (gammaI*beta)/(n+gammaI-i)
}

arg <- function(n0, n, K, k, beta, gamma, F) {
  n1 <- n-n0

  alphaP <- function(r) {
    ## sapply(r, FUN=function(ii) { 
    ##   alphaPrime(mTilde(ii, n1, gamma), beta, gamma)
    ## })
    alphaPrime(mTilde(r, n1, gamma), beta, gamma)
  }
  alphaP1 <- function(r) {
    ## sapply(r, FUN=function(ii) { 
    ##   alphaPrime(mTilde(ii, n1, gamma)-1, beta, gamma)
    ## })
    alphaPrime(mTilde(r, n1, gamma)-1, beta, gamma)
  }
  
  r1 <- alphaP(k-1)/k

  ## if (FALSE) {
  ##   r2 <- 0
  ##   rr <- k
  ##   while (rr <= K) {
  ##     deltaRR <- alphaP(rr) - alphaP(rr-1)
  ##     r2 <- r2 + deltaRR/rr
  ##     rr <- rr+1
  ##   }
  ##   rm(rr, deltaRR)
  ##   ## stopifnot(abs(r2-r2.1)<1e-10)
  ## }

  rr <- seq(from=k, length=K-k+1)
  deltaRR <- alphaP(rr) - alphaP(rr-1)
  r2 <- sum(deltaRR/rr)
  
  r3 <- 0

  rr <- seq(from=K+1, length=n0-K)
  r3 <- 0
  if (length(rr)) {
    apr <- alphaP(rr)
    ap1r <- alphaP1(rr)
    deltaRR <- apr - alphaP(rr-1)
  
    deltaSS <- sapply(seq(along=rr), FUN=function(ll) {
      r <- rr[ll]
      ss <- seq(from=r+1, length=n0-r)
      aps <- alphaP(ss)
      ap1s <- alphaP1(ss)
      GSS <- F(apr[ll], aps) - F(ap1r[ll], aps) - F(apr[ll], ap1s) + F(ap1r[ll], ap1s)  ## G(stuff)
      sum(GSS/ss)
    })

    delta2RR <- F(apr, apr) - F(apr, ap1r)
    r3 <- sum(deltaRR/rr^2 + (n0-1)*deltaSS/rr + (n0-1)*delta2RR/rr^2)
  }
  
##   if (FALSE) {
##   rr <- K+1
##   while(rr < n0) {
##     apr <- alphaP(rr)
##     ap1r <- alphaP1(rr)
##     deltaRR <- apr - alphaP(rr-1)

##     deltaSS <- 0
##     ss <- rr+1
##     while (ss < n0) {
##       aps <- alphaP(ss)
##       ap1s <- alphaP1(ss)
##       GSS <- F(apr, aps) - F(ap1r, aps) - F(apr, ap1s) + F(ap1r, ap1s)  ## G(stuff)
##       deltaSS <- deltaSS + GSS/ss
##       ss <- ss+1
##     }
    
##     delta2RR <- F(apr, apr) - F(apr, ap1r)
##     r3 <- r3 + deltaRR/rr^2 + (n0-1)*deltaSS/rr + (n0-1)*delta2RR/rr^2
##     rr <- rr+1
##   }
## }
  
  n0*(r1 + r2 + r3)
}

C3SU <- function(beta, k, gamma, F, n) {
  n0s <- seq(k, n, by=1)
  resN0 <- sapply(n0s, FUN=function(n0) {
    Ks <- seq(k, n0, by=1)
    resK <- sapply(Ks, FUN=function(KK) {
      arg(n0, n, KK, k, beta, gamma, F)
    })
    min(resK)
  })
  max(resN0)
}

betaSU <- function(alpha, k, gamma, F, n) {
  crit <- function(beta, alpha=0) C3SU(beta, k, gamma, F, n)-alpha
  uniroot(crit, c(0, 1), alpha=0.2)
}

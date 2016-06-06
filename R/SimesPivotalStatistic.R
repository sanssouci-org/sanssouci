SimesPivotalStatistic <- function(kmaxH0) {
    ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
    m <- nrow(kmaxH0)
    pval <- 1-pnorm(kmaxH0)
    skInv <- sweep(pval, 1, m/1:m, "*")
    colMins(skInv)
}

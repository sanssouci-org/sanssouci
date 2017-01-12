library("sansSouci")

typeOfSNR <- "constant"
flavor <- "equi"
sname0 <- sprintf("%s,%s", flavor, typeOfSNR)
ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
sname <- sprintf("%s,%s", sname0, ptag)

SNRs <- rev(c(0, 1, 2, 3, 4, 5))
deps <- c(0, 0.2, 0.4) ## correlation coefficient
pi0s <- c(0.8, 0.9, 0.99, 0.999)
pi0s <- c(0.8, 0.9, 0.99)
alphas <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)

m <- 1e2
B <- 1e2
nbSimu <- 1e2





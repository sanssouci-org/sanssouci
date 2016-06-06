library("sansSouci")
Rcpp <- c(TRUE, FALSE)[1]

m <- 1e3
rho <- 0
pi0 <- 0.9
B <- 1e3
#SNR <- 2
alpha <- 0.25
nbSimu <- 1e3
typeOfSNR <- c("constantSNR", "Pareto")[1]


SNRs <- c(0, 1, 2, 3, 4, 5)
#SNRs <- "Pareto(2,2,1)"
#SNRs <- c("Pareto(0,2,1)", "Pareto(1,2,1)", "Pareto(2,2,1)"); library("actuar");
rhos <- c(0, 0.2, 0.4)
pi0s <- c(0.9, 0.99, 0.999)

source("package/inst/testScripts/BNR/testStepDown.R")

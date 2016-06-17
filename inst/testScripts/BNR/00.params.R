library("sansSouci")

m <- 1e3
B <- 1e3

rho <- 0
pi0 <- 0.9
SNR <- 2

alpha <- 0.25
nbSimu <- 1e3
typeOfSNR <- c("constant", "Pareto", "linear")[3]

flavor <- c("equi", "Mein2006")[1]

SNRs <- rev(c(0, 1, 2, 3, 4, 5))
#SNRs <- "Pareto(2,2,1)"
#SNRs <- c("Pareto(0,2,1)", "Pareto(1,2,1)", "Pareto(2,2,1)"); library("actuar");
rhos <- c(0, 0.2, 0.4)
pi0s <- c(0.8, 0.9, 0.99, 0.999)

source("package/inst/testScripts/BNR/testStepDown.R")
sname0 <- sprintf("%s,%s,%s", "stepDownEqui+power", flavor, typeOfSNR)


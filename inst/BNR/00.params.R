library("sansSouci")

m <- 1e2
B <- 1e3

alpha <- 0.25
nbSimu <- 1e2

typeOfSNR <- c("constant", "Pareto", "linear")[1]
flavor <- c("equi", "Mein2006", "Toeplitz")[1]

SNRs <- rev(c(0, 1, 2, 3, 4, 5))
#SNRs <- "Pareto(2,2,1)"
#SNRs <- c("Pareto(0,2,1)", "Pareto(1,2,1)", "Pareto(2,2,1)"); library("actuar");
if (flavor == "Toeplitz") {
    deps <- c(-2, -1, -0.5, -0.2)  ## Toeplitz exponent
} else {
    deps <- c(0, 0.2, 0.4) ## correlation coefficient
}
source("package/inst/testScripts/BNR/testStepDown.R")
sname0 <- sprintf("%s,%s,%s", "stepDownEqui+power", flavor, typeOfSNR)


pi0s <- c(0.8, 0.9, 0.99, 0.999)
pi0s <- c(0.8, 0.9, 0.99)
alphas <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)


library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")

library("sansSouci")
empCov <- function(m, rho, B, pi0=1, SNR=NA, kMax=m, alpha=0.05) {
    sim <- sansSouci::simulateEqui(m, rho, B, pi0=pi0, SNR=SNR)
    sk <- sansSouci::SimesThresholdFamily(m, kMax=kMax)
    sansSouci:::empiricalCoverage(sk(alpha), sim$X0)
}
m <- 1e3
B <- 1e4
rhos <- c(0, 0.1, 0.2, 0.4, 0.8)
alpha <- 0.2

res <- listenv()
plan(remote, workers = rep("bernoulli", 20))
reps <- 100
for (rep in seq_len(reps)) {
    res[[rep]] %<-% {
        sapply(rhos, FUN=function(rho) empCov(m, rho, B, alpha=alpha))
    }
}

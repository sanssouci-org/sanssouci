## For Figures 2 and 3 in the BNR paper
library("future")
library("listenv")
library("sansSouci")
library("plyr")  ## TODO: drop this dependency

##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
#plan(remote, workers = rep("bernoulli", 100))

plan(multiprocess, workers = 3)

options("future.wait.times" = 1e4)

## setup
m <- 1e3
kMaxs <- c(1, 2, 10, m)
B <- 1e3
nbSimu <- 10

alpha <- 0.2

rhos <- seq(from = 0, to = 1, by = 0.1)

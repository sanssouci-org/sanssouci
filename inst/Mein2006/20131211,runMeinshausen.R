m <- 1000
pi0 <- 0.7

p <- 0.5

flavor <- "2006"

ns <- c(20, 40, 60, 80)
rhos <- c(0, 0.4)[1]
nbSimu <- 100

for (n in ns) {
  print(n)
  for (rho in rhos) {
    print(rho)
    source("package/inst/testScripts/Mein2006/20131211,Meinshausen,x100.R")
  }
}

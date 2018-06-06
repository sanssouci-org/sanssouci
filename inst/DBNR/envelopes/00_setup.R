ms <- 12800
ss <- c(10, 20, 50, 100)[4]
ds <- c(0.5, 0.75, 0.9, 1)
barmus <- c(2, 3, 4, 5)
K1s <- c(1, 4, 8, 16)[3]
groupeds <- c(TRUE, FALSE)
settings <- c("const", "gauss", "poisson")

configs <- expand.grid(
    m = ms,
    s = ss,
    K1 = K1s,
    d = ds,
    barmu = barmus,
    grouped = groupeds,
    setting = settings, 
    stringsAsFactors = FALSE)


repl <- 2  ## number of replications

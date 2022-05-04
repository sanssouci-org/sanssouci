library("sansSouci")
# file <- system.file("DBNR/envelopes/simu.hulk.R", package = "sansSouci", mustWork = TRUE)
# source(file)
source("inst/DBNR/envelopes/guillermo_simu.hulk.R")

ms <- 12800
ss <- rev(c(10, 50, 100, 200, 400)) # taille d'une feuille
ds <- c(0.5, 0.75, 0.9, 1) # proportion de signal dans les feuilles actives
barmus <- c(2, 3, 4, 5)[1:3]
K1s <- c(1, 4, 8, 16, 32) # nombre de feuilles avec du signal/feuilles actives
groupeds <- c(TRUE, FALSE)
settings <- c("const", "gauss", "poisson", "rgauss")[4] #settings == "rgauss"

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

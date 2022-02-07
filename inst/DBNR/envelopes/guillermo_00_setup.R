library("sansSouci")
# file <- system.file("DBNR/envelopes/simu.hulk.R", package = "sansSouci", mustWork = TRUE)
# source(file)
setwd("~/Documents/sanssouci")
source("inst/DBNR/envelopes/guillermo_simu.hulk.R")

ms <- 12800 # nombre total d'hypothèses, 12800 dans le papier
ss <- rev(c(10, 50, 100, 200, 400))[3] # taille d'une feuille, 100 dans le papier
ds <- c(0.5, 0.75, 0.9, 1)[-1] # proportion de signal dans les feuilles actives, 0.75, 0.9, 1 dans la papier
barmus <- c(2, 3, 4, 5)[1:3] # "taille de référence" du signal quand non-nul, 2, 3, 4 dans le papier
K1s <- c(1, 4, 8, 16, 32)[3] # nombre de feuilles avec du signal/feuilles actives, 8 dans le papier
groupeds <- c(TRUE, FALSE)[1] # coller ensemble ou pas les feuilles avec du signal, TRUE dans le papier
settings <- c("const", "gauss", "poisson", "rgauss")[4] # forme des espérances des stats de test, rgauss dans le papier

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

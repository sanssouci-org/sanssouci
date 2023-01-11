# set.seed(2)

source("inst/DBNR/envelopes-hybrid/hybrid_00_setup.R")
repl <- 20 # 10 dans le papier ? ## number of replications
source("inst/DBNR/envelopes-hybrid/hybrid_01_run-simulations.R")
source("inst/DBNR/envelopes-hybrid/hybrid_02_read-results.R")
source("inst/DBNR/envelopes-hybrid/hybrid_04_ggplot-paper.R")
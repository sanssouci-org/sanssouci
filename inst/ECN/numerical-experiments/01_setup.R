library("sanssouci")
library("GSEABenchmarkeR")
library("future.apply")
library("cherry")
library("tidyr")
library("dplyr")
library("ggplot2")

# set.seed(0xBEEF) # for reproducibility

alphas <- seq(from = 0, to = 1, by = 0.05)  # target JER level
B <- 1e3          # number of permutations for adaptive methods
nb_exp <- 1e3     # number of experiments

# future::availableCores() to know available 'workers'
plan(multisession, workers = 40) 

technology <- c("microarray", "RNAseq")[1]

# create data set and experiment parameters
if (technology == "microarray") {
    rowTestFUN <- sanssouci::rowWelchTests
    
    geo2kegg <-  R.cache::memoizedCall(loadEData,"geo2kegg")
    ds_name <- "GSE19188"
    rawData <- R.cache::memoizedCall(maPreproc, geo2kegg[ds_name])[[1]]
    X <- SummarizedExperiment::assays(rawData)$exprs
    cats <- SummarizedExperiment::colData(rawData)
    ww <- match(cats$Sample, base::colnames(X))
    groups <- cats$GROUP[ww]
    
    pi0s <- c(0.8, 0.95, 1)
    SNR <- c(0, 1, 2, 5)
    SNR_FUN <- "+"
    probs <- c(0.5, 0.25, 0.1)
    
    X0 <- X[, groups == 0]
    
} else if (technology == "RNAseq") {
    rowTestFUN <- sanssouci::rowWilcoxonTests
    
    data("RNAseq_blca", package = "sanssouci.data")
    ds_name <- "BLCA"
    X <- RNAseq_blca
    groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
    rm(RNAseq_blca)
    
    pi0s <- c(0.8, 0.95, 1)
    SNR <- c(1, 2, 3)
    SNR_FUN <- "*"
    probs <- 0.5
    
    X0 <- X[, groups == 0]
    
    # filter out unexpressed genes
    if (technology == "RNAseq"){
        BLCA0 <- X0/colSums(X)*1e6
        ww <- which(rowQuantiles(BLCA0, prob = 0.75) < 5)
        if (length(ww) != 0){
            X0 <- X0[-ww, ]
        }
    }
    
}



table(groups)

dim(X0)
m <- nrow(X0)

configs <- expand.grid(SNR = SNR,
                       pi0 = pi0s,
                       SNR_FUN = SNR_FUN,
                       prob = probs)
seq_configs <- 1:nrow(configs)

path <- sprintf("results/diff-expr_%s", technology)
dir.create(path, showWarnings = FALSE, recursive = TRUE)


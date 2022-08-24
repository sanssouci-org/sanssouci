library("ggplot2")
library("tidyr")
library("sanssouci")
library("matrixStats")

data("RNAseq_blca", package = "sanssouci.data")
Y <- RNAseq_blca
groups <- ifelse(colnames(RNAseq_blca) == "III", 1, 0)
rm(RNAseq_blca)
dim(Y)

CPM <- Y/colSums(Y)*1e6
ww <- which(rowQuantiles(CPM, prob = 0.75) < 5)
if (length(ww) != 0){
    Y <- Y[-ww, ]
}
rm(CPM, ww)
dim(Y)

set.seed(19012001)

alpha <- 0.1
obj <- SansSouci(Y = log(1 + Y), groups = groups)
res <- fit(obj, B = 1000, alpha = alpha, family = "Simes", 
           rowTestFUN = rowWilcoxonTests)

res_singlestep <- fit(obj, B = 1000, alpha = alpha, family = "Simes", 
                      rowTestFUN = rowWilcoxonTests, max_steps_down = 0)

res_Simes <- fit(obj, B = 0, family = "Simes", alpha = alpha, 
                 rowTestFUN = rowWilcoxonTests) ## B=0 => no calibration!

m <- dim(Y)[1]
d <- hommel::discoveries(hommel::hommel(pValues(res)))
pi0_hat <- 1 - d/m

res_ARI <- fit(obj, B = 0, family = "Simes", alpha = alpha/pi0_hat, 
               rowTestFUN = rowWilcoxonTests)


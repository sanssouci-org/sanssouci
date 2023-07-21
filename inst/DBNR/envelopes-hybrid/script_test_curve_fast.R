library(sanssouci)
library(microbenchmark)

set.seed(12)

# C <- list(
# 	list(c(2, 5), c(8, 15)), 
# 	list(c(3, 5), c(8, 10), c(12, 15)), 
# 	list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15)),
# 	list(c(8, 8), c(9, 9), c(13, 13), c(14, 15))
# )
# ZL <- list(
# 	c(4, 8), 
# 	c(3, 3, 4), 
# 	c(2, 2, 1, 1, 2),
# 	c(1, 1, 1, 2)
# )
# leaf_list <- as.list(1:16)
# completed <- forest.completion(C, ZL, leaf_list)

m <- 16
example <- dyadic.from.height(m, 3, 2)
leaf_list <- example$leaf_list
C <- example$C

pval_idx <- unlist(leaf_list)
pval <- runif(n = max(pval_idx))

alpha = 0.05
method <- zeta.HB

pval[1:4] <- 1e-15

print("The pvalues are:")
print(pval)

zetas <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
print("The zetas are:")
print(zetas)
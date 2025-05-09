library(sanssouci)

set.seed(12)

m <- 16
example <- dyadic.from.height(m, 3, 2)
leaf_list <- example$leaf_list
C <- example$C

# pval_idx <- unlist(leaf_list)
# pval <- runif(n = max(pval_idx))
pval <- runif(m)

alpha = 0.05
method <- zeta.HB

pval[1:4]<-1e-15
pval[5] <- 0.002

print("The pvalues are:")
print(pval)

print("Without refining:")
zetas_no <- zetas.tree(C, leaf_list, method, pval, alpha, refine = FALSE, verbose = TRUE)
print("The zetas are:")
print(zetas_no)
print("With refining:")
zetas_yes <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = TRUE)
print("The zetas are:")
print(zetas_yes)

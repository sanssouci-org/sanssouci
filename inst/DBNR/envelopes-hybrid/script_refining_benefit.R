set.seed(12)

example <- dyadic.from.height(16,3,2)
leaf_list <- example$leaf_list
C <- example$C

pval_idx <- unlist(leaf_list)
pval <- runif(n=max(pval_idx))

alpha=0.05
method<-zeta.HB

pval[1:4]<-1e-15

print("The pvalues are:")
print(pval)

print("Without refining:")
zetas_no <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine=FALSE, verbose=TRUE)
print("The zetas are:")
print(zetas_no)
print("With refining:")
zetas_yes <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine=TRUE, verbose=TRUE)
print("The zetas are:")
print(zetas_yes)

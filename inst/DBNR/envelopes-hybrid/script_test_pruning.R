# TODO BEFORE MERGE: rename zetas.tree
library(sanssouci)
library(microbenchmark)

set.seed(12)

m <- 16
example <- dyadic.from.height(m, 3, 2)
leaf_list <- example$leaf_list
C <- example$C

pval_idx <- unlist(leaf_list)
pval <- runif(n = max(pval_idx))

alpha <- 0.05
method <- zeta.HB

pval[1:4] <- 1e-15

print("The pvalues are:")
print(pval)

zetas <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
print("The zetas are:")
print(zetas)

pruned <- pruning(C, zetas, leaf_list, TRUE)
print(paste0("V^*({1, ..., ", m, "})=", pruned$VstarNm))
print("New structure:")
print(pruned$C)
print("New associated zetas:")
print(pruned$ZL)

print("Testing the pruned structure")
print(paste0("V^*({1, ..., ", m, "})=", V.star(1:m, pruned$C, pruned$ZL, leaf_list)))
print(paste0("V^*({1, ..., ", 5, "})=", V.star(1:5, pruned$C, pruned$ZL, leaf_list)))
print(paste0("V^*({1, ..., ", 8, "})=", V.star(1:8, pruned$C, pruned$ZL, leaf_list)))
print(paste0("V^*({9, ..., ", 16, "})=", V.star(9:16, pruned$C, pruned$ZL, leaf_list)))


print("Changing for an arbitrary zeta tree to run more tests")
arbitrary_zetas <- list(10,
							c(4, 6),
							c(0, 4, 3 , 4)
							)
print("The zetas now are:")
print(arbitrary_zetas)

pruned2 <- pruning(C, arbitrary_zetas, leaf_list, TRUE)
print(paste0("V^*({1, ..., ", m, "})=", pruned2$VstarNm))
print("New structure:")
print(pruned2$C)
print("New associated zetas:")
print(pruned2$ZL)

print("Testing the pruned structure")
print(paste0("V^*({1, ..., ", m, "})=", V.star(1:m, pruned2$C, pruned2$ZL, leaf_list)))
print(paste0("V^*({1, ..., ", 5, "})=", V.star(1:5, pruned2$C, pruned2$ZL, leaf_list)))
print(paste0("V^*({1, ..., ", 8, "})=", V.star(1:8, pruned2$C, pruned2$ZL, leaf_list)))
print(paste0("V^*({9, ..., ", 16, "})=", V.star(9:16, pruned2$C, pruned2$ZL, leaf_list)))

print("Now we test the computation time enhancement provided by the pruning, on S={1, ..., m} and S={5}")
pruned3 <- pruning(C, zetas, leaf_list, FALSE)

mbench <- microbenchmark(no_pruning = V.star(1:m, C, zetas, leaf_list), 
							 pruning = V.star(1:m, pruned3$C, pruned3$ZL, leaf_list),
							 super_pruning = V.star(1:m, pruned$C, pruned$ZL, leaf_list),
							 check = "equal",
							 times = 1000
							 )
print(mbench)

mbench <- microbenchmark(no_pruning = V.star(5, C, zetas, leaf_list), 
												 pruning = V.star(5, pruned3$C, pruned3$ZL, leaf_list),
												 super_pruning = V.star(5, pruned$C, pruned$ZL, leaf_list),
												 check = "equal",
												 times = 1000
)
print(mbench)




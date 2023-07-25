# TODO BEFORE MERGE: rename zetas.tree
library(sanssouci)
library(microbenchmark)

set.seed(12)

#########################################################
################ FIRST EXAMPLE
#########################################################
m <- 20
C <- list(
	list(c(2, 5), c(8, 15), c(16, 19)),
	list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
	list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15), c(17, 17), c(18, 19)),
	list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
)
ZL <- list(
	c(4, 8, 4),
	c(3, 3, 4, 1, 3),
	c(2, 2, 1, 1, 2, 1, 2),
	c(1, 1, 1, 2, 1, 1)
)
leaf_list <- as.list(1:m)
completed <- forest.completion(C, ZL, leaf_list)
C <- completed$C
ZL <- completed$ZL

# print("Printing assets")
# print(C)
# print(ZL)

pruned <- pruning(C, ZL, leaf_list, super.prune = FALSE)

# print("Printing pruned assets")
# print(pruned$C)
# print(pruned$ZL)

print("test delete.gaps")
gaps.deleted <- delete.gaps(pruned$C, pruned$ZL, leaf_list)
print(gaps.deleted)


# print("Testing that the three versions of curve.V.star have the same output:")
# perm <- 1:m
# print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
# print(curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, pruning = FALSE, is.pruned = TRUE))
# print(curve.V.star.forest.fast(perm, gaps.deleted$C, gaps.deleted$ZL, leaf_list, pruning = FALSE, is.pruned = TRUE))


print("Comparing execution times:")
mbench <- microbenchmark(fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list, is.complete = TRUE),
												 fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
												 fast.pruned.no.gaps =  curve.V.star.forest.fast(perm, gaps.deleted$C, gaps.deleted$ZL, leaf_list, is.pruned = TRUE),
												 times=100, check="equal")
print(mbench)

#########################################################
################ SECOND EXAMPLE
#########################################################
pow <- 10
m <- 2 ^ pow
example <- dyadic.from.height(m, pow, 2)
leaf_list <- example$leaf_list
C <- example$C
pval <- runif(n = m)
alpha <- 0.05
method <- zeta.trivial
ZL <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
completed <- forest.completion(C, ZL, leaf_list)
C <- completed$C
ZL <- completed$ZL

# print("Printing assets")
# print(C)
# print(ZL)

pruned <- pruning(C, ZL, leaf_list, super.prune = FALSE)

# print("Printing pruned assets")
# print(pruned$C)
# print(pruned$ZL)

print("test delete.gaps")
gaps.deleted <- delete.gaps(pruned$C, pruned$ZL, leaf_list)
print(gaps.deleted$ZL)


# print("Testing that the three versions of curve.V.star have the same output:")
# perm <- 1:m
# print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
# print(curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, pruning = FALSE, is.pruned = TRUE))
# print(curve.V.star.forest.fast(perm, gaps.deleted$C, gaps.deleted$ZL, leaf_list, pruning = FALSE, is.pruned = TRUE))


print("Comparing execution times:")
mbench <- microbenchmark(fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list, is.complete = TRUE),
												 fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
												 fast.pruned.no.gaps =  curve.V.star.forest.fast(perm, gaps.deleted$C, gaps.deleted$ZL, leaf_list, is.pruned = TRUE),
												 times=100, check="equal")
print(mbench)



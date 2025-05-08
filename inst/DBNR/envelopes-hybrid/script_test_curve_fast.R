library(sanssouci)
library(microbenchmark)

set.seed(12)

#########################################################
################ FIRST EXAMPLE
#########################################################
m <- 16
example <- dyadic.from.height(m, 3, 2)
leaf_list <- example$leaf_list
C <- example$C
pval <- runif(n = m)
alpha <- 0.05
method <- zeta.HB
pval[13:16] <- 1e-15
ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)

print("The pvalues are:")
print(pval)
print("The zetas are:")
print(ZL)

print("Testing that the four versions of curve.V.star have the same output:")
perm <- 1:m
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = TRUE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = TRUE))

super.pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE)
pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE)

print("Comparing execution times:")
mbench <- microbenchmark(naive.not.pruned = curve.V.star.forest.naive(perm, C, ZL, leaf_list),
												 naive.pruned = curve.V.star.forest.naive(perm, super.pruned$C, super.pruned$ZL, leaf_list),
												 fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list),
												 fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
												 times=10, check="equal")
print(mbench)

print("Testing that the four versions of curve.V.star have the same output:")
perm <- m:1
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = TRUE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = TRUE))

super.pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE)
pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE)

print("Comparing execution times:")
mbench <- microbenchmark(naive.not.pruned = curve.V.star.forest.naive(perm, C, ZL, leaf_list),
												 naive.pruned = curve.V.star.forest.naive(perm, super.pruned$C, super.pruned$ZL, leaf_list),
												 fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list),
												 fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
												 times=10, check="equal")
print(mbench)


#########################################################
################ SECOND EXAMPLE
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

print("The zetas are:")
print(ZL)

print("Testing that the four versions of curve.V.star have the same output:")
perm <- 1:m
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = TRUE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = TRUE))

super.pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE)
pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE)

print("Comparing execution times:")
mbench <- microbenchmark(naive.not.pruned = curve.V.star.forest.naive(perm, C, ZL, leaf_list),
												 naive.pruned = curve.V.star.forest.naive(perm, super.pruned$C, super.pruned$ZL, leaf_list),
												 fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list),
												 fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
												 times=10, check="equal")
print(mbench)


#########################################################
################ THIRD EXAMPLE
#########################################################
pow <- 10
m <- 2 ^ pow
example <- dyadic.from.height(m, pow, 2)
leaf_list <- example$leaf_list
C <- example$C
pval <- runif(n = m)
alpha <- 0.05
method <- zeta.trivial
ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)

print("The pvalues are:")
print(pval)
print("The zetas are:")
print(ZL)

print("Testing that the four versions of curve.V.star have the same output:")
perm <- 1:m
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.naive(perm, C, ZL, leaf_list, pruning = TRUE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = FALSE))
print(curve.V.star.forest.fast(perm, C, ZL, leaf_list, pruning = TRUE))

super.pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE)
pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE)

print("Comparing execution times:")
mbench <- microbenchmark(naive.not.pruned = curve.V.star.forest.naive(perm, C, ZL, leaf_list),
                         naive.pruned = curve.V.star.forest.naive(perm, super.pruned$C, super.pruned$ZL, leaf_list),
                         fast.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list),
                         fast.pruned = curve.V.star.forest.fast(perm, pruned$C, pruned$ZL, leaf_list, is.pruned = TRUE),
                         times=10, check="equal")
print(mbench)

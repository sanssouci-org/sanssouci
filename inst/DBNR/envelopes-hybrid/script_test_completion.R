library(sanssouci)

C <- list(
	list(c(2, 5), c(8, 15)), 
	list(c(3, 5), c(8, 10), c(12, 15)), 
	list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15)),
	list(c(8, 8), c(9, 9), c(13, 13), c(14, 15))
)
ZL <- list(
	c(4, 8), 
	c(3, 3, 4), 
	c(2, 2, 1, 1, 2),
	c(1, 1, 1, 2)
)
leaf_list <- as.list(1:16)

desired_C_outcome <- list(
	list(c(1, 1), c(2, 5), c(6, 6), c(7, 7), c(8, 15), c(16, 16)), 
	list(c(2, 2), c(3, 5), c(8, 10), c(11, 11), c(12, 15)), 
	list(c(3, 3), c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15)),
	list(c(4, 4), c(5, 5), c(8, 8), c(9, 9), c(13, 13), c(14, 15)),
	list(c(14, 14), c(15, 15))
)
desired_ZL_outcome <- list(
	c(1, 4, 1, 1, 8, 1), 
	c(1, 3, 3, 1, 4), 
	c(1, 2, 2, 1, 1, 2),
	c(1, 1, 1, 1, 1, 2),
	c(1, 1)
)

completed <- forest.completion(C, ZL, leaf_list)

print(completed$C)
print(completed$ZL)

print(all.equal(completed$C, desired_C_outcome) && all.equal(completed$ZL, desired_ZL_outcome))

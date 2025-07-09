#! /usr/bin/env Rscript
library(sanssouci)
library(microbenchmark)
library(r2r)

set.seed(12)


pow <- 10
alpha <- 0.05

vec_factor <- c(1, 2, 4, 8)
record_time_tree <- rep(0, length(vec_factor))
record_time_part <- rep(0, length(vec_factor))
method <- zeta.trivial
n_repl <- 100

for (i in 1:length(vec_factor)){

  factor <- vec_factor[i]
  m <- (2 ^ pow) * factor
  example <- dyadic.from.height(m, pow - 1, 2)
  leaf_list <- example$leaf_list
  C <- example$C
  H <- length(C)
  pval <- 1 - pnorm(rnorm(n = m))
  ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)

  C2 <- list(C[[H]])
  ZL2 <- list(ZL[[H]])

  print(paste0("m = ", m))
  print(paste0("K pour part = ", nb.elements(C)))
  print(paste0("K pour tree = ", nb.elements(C2)))
  print(paste0("N = ", length(leaf_list)))

  print("Comparing execution times:")
  mbench <- microbenchmark(tree = V.star2(1:m, C, ZL, leaf_list),
                           part = V.star2(1:m, C2, ZL2, leaf_list),
                           times = n_repl, check="equal")
  print(mbench, unit="milliseconds")
  record_time_tree[i] <- summary(mbench, unit="milliseconds")$median[1]
  record_time_part[i] <- summary(mbench, unit="milliseconds")$median[2]
}
plot((2 ^ pow) * vec_factor, record_time_tree, ylim = c(0, 7))
points((2 ^ pow) * vec_factor, record_time_part, col="red")

vec_factor <- c(1, 2, 4, 8, 16)
# factor <- 1
stacks <- 10
record_time_part <- rep(0, length(vec_factor))
record_time_part_stacked <- rep(0, length(vec_factor))
record_time_part_half <- rep(0, length(vec_factor))
record_time_part_stacked_half <- rep(0, length(vec_factor))
method <- zeta.trivial
n_repl <- 100

for (i in 1:length(vec_factor)){
  
  factor <- vec_factor[i]
  m <- (2 ^ pow) * factor
  example <- dyadic.from.height(m, pow - 1, 2)
  leaf_list <- example$leaf_list
  C <- example$C
  H <- length(C)
  pval <- 1 - pnorm(rnorm(n = m))
  ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
  
  C2 <- list(C[[H]])
  ZL2 <- list(ZL[[H]])
  stackC <- vector("list", length = stacks)
  stackZL <- vector("list", length = stacks)
  for (j in 1:stacks){
    stackC[[j]] = C2[[1]]
    stackZL[[j]] = ZL2[[1]]
  }
  
  print(paste0("m = ", m))
  print(paste0("K pour tree = ", nb.elements(C)))
  print(paste0("K pour part = ", nb.elements(C2)))
  print(paste0("K pour part stacked = ", nb.elements(stackC)))
  print(paste0("N = ", length(leaf_list)))
  
  print("Comparing execution times:")
  mbench <- microbenchmark(#tree = V.star2(1:m, C, ZL, leaf_list),
                           part = V.star2(1:m, C2, ZL2, leaf_list),
                           part_stacked = V.star2(1:m, stackC, stackZL, leaf_list),
                           times = n_repl, check="equal")
  print(mbench, unit="milliseconds")
  record_time_part[i] <- summary(mbench, unit="milliseconds")$median[1]
  record_time_part_stacked[i] <- summary(mbench, unit="milliseconds")$median[2]
}
for (i in 1:length(vec_factor)){
  
  factor <- vec_factor[i]
  m <- (2 ^ pow) * factor
  example <- dyadic.from.height(m, pow - 2, 2)
  leaf_list <- example$leaf_list
  C <- example$C
  H <- length(C)
  pval <- 1 - pnorm(rnorm(n = m))
  ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
  
  C2 <- list(C[[H]])
  ZL2 <- list(ZL[[H]])
  stackC <- vector("list", length = stacks)
  stackZL <- vector("list", length = stacks)
  for (j in 1:stacks){
    stackC[[j]] = C2[[1]]
    stackZL[[j]] = ZL2[[1]]
  }
  
  print(paste0("m = ", m))
  print(paste0("K pour tree = ", nb.elements(C)))
  print(paste0("K pour part = ", nb.elements(C2)))
  print(paste0("K pour part stacked = ", nb.elements(stackC)))
  print(paste0("N = ", length(leaf_list)))
  
  print("Comparing execution times:")
  mbench <- microbenchmark(#tree = V.star2(1:m, C, ZL, leaf_list),
    part = V.star2(1:m, C2, ZL2, leaf_list),
    part_stacked = V.star2(1:m, stackC, stackZL, leaf_list),
    times = n_repl, check="equal")
  print(mbench, unit="milliseconds")
  record_time_part_half[i] <- summary(mbench, unit="milliseconds")$median[1]
  record_time_part_stacked_half[i] <- summary(mbench, unit="milliseconds")$median[2]
}
plot((2 ^ pow) * vec_factor, record_time_part, ylim = c(0, 12))
points((2 ^ pow) * vec_factor, record_time_part_stacked, col="red")
points((2 ^ pow) * vec_factor, record_time_part_half, col="blue")
points((2 ^ pow) * vec_factor, record_time_part_stacked_half, col="green")


# factor <- vec_factor[1]
# m <- (2 ^ pow) * factor
# example <- dyadic.from.height(m, pow - 1, 2)
# leaf_list <- example$leaf_list
# nb_leaves <- length(leaf_list)
# 
# C <- example$C
# pval <- 1 - pnorm(rnorm(n = m))
# ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
# S <- 1:m
# 
# print(paste0("m = ", m))
# print(paste0("K = ", nb.elements(C)))
# print(paste0("N = ", length(leaf_list)))
# 
# print("Comparing execution times:")
# mbench <- microbenchmark(V.star = V.star(S, C, ZL, leaf_list),
#                          V.star2 = V.star2(S, C, ZL, leaf_list),
#                          times = n_repl, check="equal")
# print(mbench)



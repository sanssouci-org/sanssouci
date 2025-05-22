#! /usr/bin/env Rscript
library(sanssouci)
library(microbenchmark)
library(r2r)

set.seed(12)

V.star2 <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- sum(leaf_list[[i]][length(leaf_list[[i]])] >= S)
  }
  Vec <- Vec - c(0, Vec[1:(nb_leaves-1)])
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  # this initialization also takes care of the minima
  # between \zeta_k and card(S inter R_k)
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions > 0) {
      for (k in 1:nb_regions) {
        Rk <- C[[h]][[k]]
        sum_succ <- sum(Vec[Rk[1]:Rk[2]])
        res <- min(ZL[[h]][k], sum_succ)
        Vec[Rk[1]:Rk[2]] <- 0
        Vec[Rk[1]] <- res
      }
    }
  }
  return(sum(Vec))
}

pow <- 10
alpha <- 0.05


vec_factor <- c(1, 2, 4, 8)
record_time <- rep(0, length(vec_factor))
record_time_init <- rep(0, length(vec_factor))
method <- zeta.trivial

n_repl <- 100


for (i in 1:length(vec_factor)){
  
  factor <- vec_factor[i]
  # factor <- 1
  m <- (2 ^ pow) * factor
  example <- dyadic.from.height(m, pow - 1, 2)
  leaf_list <- example$leaf_list
  C <- example$C
  H <- length(C)
  pval <- 1 - pnorm(rnorm(n = m))
  ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
  
  print(paste0("m = ", m))
  print(paste0("K pour part = ", nb.elements(C)))
  print(paste0("K pour tree = ", nb.elements(C2)))
  print(paste0("N = ", length(leaf_list)))
  
  
  C2 <- list(C[[H]])
  ZL2 <- list(ZL[[H]])

  print("Comparing execution times:")
  mbench <- microbenchmark(tree = V.star2(1:m, C, ZL, leaf_list),
                           part = V.star2(1:m, C2, ZL2, leaf_list),
                           times = n_repl, check="equal")
  print(mbench, unit="milliseconds")
  record_time[i] <- summary(mbench, unit="milliseconds")$median[1]
  record_time_init[i] <- summary(mbench, unit="milliseconds")$median[2]
}




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



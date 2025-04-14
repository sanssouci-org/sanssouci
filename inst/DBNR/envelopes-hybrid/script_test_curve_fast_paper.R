#! /usr/bin/env Rscript
library(sanssouci)
library(microbenchmark)

source("inst/DBNR/envelopes-hybrid/new_impl_fast.R")

set.seed(12)

get_groups <- function(list_groups, leaf_list){
  res <- list()
  for (group in list_groups){
    res <- c(res, list(leaf_list[[group[1]]][1]:rev(leaf_list[[group[2]]])[1]))
  }
  return(list(res = res, total = unlist(res)))
}

pow <- 10
alpha <- 0.05

vec_factor <- c(1, 10)
vec_n_repl <- c(100, 10)
vec_method <- c(zeta.DKWM, zeta.trivial)

for (i in 1:4){
  
  factor <- vec_factor[floor((i+1)/2)]
  n_repl <- vec_n_repl[floor((i+1)/2)]
  method <- vec_method[[i %% 2 + 1]]
  
  m <- (2 ^ pow) * factor
  example <- dyadic.from.height(m, pow, 2)
  leaf_list <- example$leaf_list
  C <- example$C
  signal <- 4
  mu <- rep(0, m)
  H1 <- get_groups(example$C[[6]][c(1, 5, 9, 10)], leaf_list)$total
  mu[H1] <- signal
  
  pval <- 1 - pnorm(mu + rnorm(n = m))
  
  ZL <- zetas.tree(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)
  
  pruned.no.gaps <- pruning(C, ZL, leaf_list, prune.leafs = FALSE, delete.gaps = TRUE)
  
  print(m)
  print(nb.elements(C))
  print(nb.elements(pruned.no.gaps$C))
  
  perm <- 1:m
  
  print("Comparing execution times:")
  mbench <- microbenchmark(#naive.not.pruned = curve.V.star.forest.naive(perm, C, ZL, leaf_list),
                           #naive.pruned = curve.V.star.forest.naive(perm, pruned.no.gaps$C, pruned.no.gaps$ZL, leaf_list),
                           fast13.not.pruned = curve.V.star.forest.fast(perm, C, ZL, leaf_list),
                           fast13.pruned = curve.V.star.forest.fast(perm, pruned.no.gaps$C, pruned.no.gaps$ZL, leaf_list, is.pruned = TRUE),
                           fast14hypcol.not.pruned = curve.V.star.forest.fast.14hypcol(perm, C, ZL, leaf_list),
                           fast14hypcol.pruned = curve.V.star.forest.fast.14hypcol(perm, pruned.no.gaps$C, pruned.no.gaps$ZL, leaf_list, is.pruned = TRUE),
                           fast14hyprow.not.pruned = curve.V.star.forest.fast.14hyprow(perm, C, ZL, leaf_list),
                           fast14hyprow.pruned = curve.V.star.forest.fast.14hyprow(perm, pruned.no.gaps$C, pruned.no.gaps$ZL, leaf_list, is.pruned = TRUE),
                           times = n_repl, check = "equal")
  print(mbench)
  #write.csv(mbench, paste0("benchmark_0", i, ".csv"), row.names = F)
}

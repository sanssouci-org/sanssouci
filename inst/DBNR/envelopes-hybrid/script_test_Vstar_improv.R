# TODO BEFORE MERGE: rename zetas.tree
library(sanssouci)
library(microbenchmark)

V.star.1 <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- sum(S %in% leaf_list[[i]])
  }
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions>0) {
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        if (Chj[1]==Chj[2]) { # means this is an atom, no need to compute 
          # len_inter given that we already did it during initialization,
          # furthermore there are no successors
          len_inter <- Vec[Chj[1]]
          res <- min(ZL[[h]][j], len_inter)
        } else {
          region_vector <- unlist(leaf_list[Chj[1]:Chj[2]])
          len_inter <- sum(S %in% region_vector)
          sum_succ <- sum(Vec[Chj[1]:Chj[2]]) 
          res <- min(ZL[[h]][j], len_inter, sum_succ)
        }
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(sum(Vec))
}

V.star.2 <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- sum(S %in% leaf_list[[i]])
  }
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions>0) {
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        
        #region_vector <- unlist(leaf_list[Chj[1]:Chj[2]])
        #len_inter <- sum(S %in% region_vector)
        sum_succ <- sum(Vec[Chj[1]:Chj[2]]) 
        # res <- min(ZL[[h]][j], len_inter, sum_succ)
        res <- min(ZL[[h]][j], sum_succ)
        
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(sum(Vec))
}


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


micr <- microbenchmark(
  V1 = V.star.1(1:m, C, ZL, leaf_list),
  V2 = V.star.2(1:m, C, ZL, leaf_list),
  check = "equal"
)

print(micr)

#########################################################
################ SECOND EXAMPLE
#########################################################

m <- 20
C <- list(
  list(c(1,2)),
  list(c(1,1), c(2,2))
)
ZL <- list(
  c(8),
  c(2, 5)
)
leaf_list <- list(1:10, 11:20)

micr <- microbenchmark(
  V1 = V.star.1(1:m, C, ZL, leaf_list),
  V2 = V.star.2(1:m, C, ZL, leaf_list),
  check = "equal"
)

print(micr)

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
ZL <- zetas.tree.no.extension(C, leaf_list, method, pval, alpha, refine = TRUE, verbose = FALSE)


micr <- microbenchmark(
  V1 = V.star.1(1:m, C, ZL, leaf_list),
  V2 = V.star.2(1:m, C, ZL, leaf_list),
  check = "equal"
)

print(micr)

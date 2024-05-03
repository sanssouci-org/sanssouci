# TODO dyadic.from.max.size doesn't exist !! and not in the examples
# dyadic.from.leaf_list is not documented and not in the examples, stranges last line in examples
#' Create a complete dyadic tree structure
#' 
#' @name dyadic
#' @details \describe{
#' \item{\code{dyadic.from.window.size}}{Dyadic tree structure from window size: the number of elements in each leaf is set to \code{s}}
#' \item{\code{dyadic.from.height}}{Dyadic tree structure from height: the total height of tree is set to \code{H}}
#' \item{\code{dyadic.from.max.size}}{Dyadic tree structure from maximum achievable height: the total height of tree is set to the maximum possible height, ie \code{floor(2 + log2(m - 1))}}
#' }
#' @examples
#' m <- 6
#' dd <- dyadic.from.window.size(m, s = 2, method = 2)
#' str(dd)
#' 
#' dd <- dyadic.from.height(m, H = 3, method = 2)
#' str(dd)
#' 
#' dd <- dyadic.from.height(m, method=2)
#' str(dd)
#'
#' leaf_list <- dd$leaf_list
NULL

#' @param leaf_list A list of leaves
#' @param method A numeric value. If \code{method == 1}, start from the leaves
#'   and group nodes of a same height 2 by 2 as long as possible. If
#'   \code{method==2}, start from the root and divide nodes in 2 nodes of equal
#'   size as long as possible
#' @return A list of lists containing the dyadic structure
#' @export
#' @rdname dyadic
dyadic.from.leaf_list <- function(leaf_list, method) {
  leaves <- length(leaf_list)
  if (method == 1) {
    inter <- seq_len(leaves)
    Ch <- mapply(c, inter, inter, SIMPLIFY = FALSE)
    C <- list(Ch)
    repeat {
      len <- length(Ch)
      if (len == 1) 
        break
      new_Ch <- list()
      j <- 1
      while (j <= len) {
        if (j == len) {
          new_Ch <- c(new_Ch, Ch[j])
        } else {
          new_Ch <- c(new_Ch, list(c(Ch[[j]][1], Ch[[j + 1]][2])))
        }
        j <- j + 2
      }
      Ch <- new_Ch
      C <- c(list(Ch), C)
    }
  } else if (method == 2) {
    Ch <- list(c(1, leaves))
    C <- list(Ch)
    continue <- TRUE
    while (continue) {
      continue <- FALSE
      oldCh <- Ch
      Ch <- list()
      len <- length(oldCh)
      for (i in seq_len(len)) {
        oi <- oldCh[[i]]
        leaves_in_node <- oi[2] - oi[1] + 1
        if (leaves_in_node > 1) {
          cut2 <- ceiling(leaves_in_node/2)
          Ch <- c(Ch, 
                  list(c(oi[1], oi[1] + cut2 - 1)), 
                  list(c(oi[1] + cut2, oi[2])))
          if (cut2 > 1) 
            continue <- TRUE
        } else {
          Ch <- c(Ch, oldCh[i])
        }
      }
      C <- c(C, list(Ch))
    }
  }
  return(C)
}

#' @param m An integer value, the number of elements in the structure
#' @param s An integer value, the number of elements in each leaf
#' @return A list with two elements:\describe{
#' \item{\code{leaf_list}}{A list of leaves}
#' \item{\code{C}}{A list of lists containing the dyadic structure}
#' }
#' @export
#' @rdname dyadic
dyadic.from.window.size <- function(m, s, method) {
  leaves <- floor(m/s)
  leaf_list <- list()
  for (l in 1:leaves) {
    leaf <- seq_len(s) + (l - 1) * s
    leaf_list <- c(leaf_list, list(leaf))
  }
  if (s * leaves < m) {
    leaf <- seq(1 + l * s, m)
    leaf_list <- c(leaf_list, list(leaf))
    # leaves<-leaves+1
  }
  C <- dyadic.from.leaf_list(leaf_list, method)
  return(list(leaf_list = leaf_list, C = C))
}

#' @param H An integer value, the desired maximal height of the tree
#' @export
#' @rdname dyadic
dyadic.from.height <- function(m, H = NULL, method) {
  if (is.null(H)) {
    H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
  }
  if (m <= 2^(H - 2)) {
    oldH <- H
    H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
    warning("H=", oldH, " is too large for m=", m, ", \nH=", oldH, " reduced to H=", H)
  }
  if (m < 2^(H - 1)) {
    leaf_list <- as.list(1:m)
  } else {
    leaf_list <- list()
    base <- m %/% 2^(H - 1)
    plus1 <- m %% 2^(H - 1)
    for (i in seq_len(plus1)) {
      leaf_list <- c(leaf_list, list(seq_len(base + 1) + (i - 1) * (base + 1)))
    }
    for (i in seq_len(2^(H - 1) - plus1)) {
      leaf_list <- c(leaf_list, list(plus1 * (base + 1) + seq_len(base) + (i - 1) * base))
    }
  }
  C <- dyadic.from.leaf_list(leaf_list, method)
  return(list(leaf_list = leaf_list, C = C))
}

#' Estimate the number of true null hypotheses among a set of p-values
#' 
#' @description
#' An upper bound on the number of true null hypotheses in the region associated to 
#' the \eqn{p}-values \code{pval}
#' is computed with confidence \code{1 - lambda}. 
#' The functions described here can be used as the \code{method} argument 
#' of [zetas.tree()].
#'
#' @param pval A vector of \eqn{p}-values
#' @param lambda A numeric value in \eqn{[0,1]}, the target error level of the test 
#' @name zeta
#' @return The number of true nulls is over-estimated as follows:
#' \describe{
#' \item{\code{zeta.DKWM}}{Inversion of the Dvoretzky-Kiefer-Wolfowitz-Massart inequality (related to the Storey estimator of the proportion of true nulls) with parameter \code{lambda}}
#' \item{\code{zeta.HB}}{Number of conserved hypotheses of the Holm-Bonferroni procedure with parameter \code{lambda}}
#' \item{\code{zeta.trivial}}{The size of the p-value set which is the trivial upper bound (\eqn{lambda} is not used)}
#' }
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. The Annals of Mathematical Statistics, pages 642-669.
#' @references Holm, S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6 (1979), pp. 65-70.
#' @references Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. The Annals of Probability, pages 1269-1283.
#' @references Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(3):479-498.
#' @examples
#' x <- rnorm(100, mean = c(rep(c(0, 2), each = 50)))
#' pval <- 1 - pnorm(x)
#' lambda <- 0.05
#' zeta.trivial(pval, lambda)
#' 
#' zeta.HB(pval, lambda)
#' 
#' zeta.DKWM(pval, lambda)
NULL

#' @rdname zeta
#' @export
zeta.HB <- function(pval, lambda) {
  m <- length(pval)
  sorted.pval <- sort(pval)
  
  thresholds <- lambda / (m - 1:m + 1)
  v <- sorted.pval - thresholds
  indexes <- which(v > 0)
  if (!length(indexes)){
    return(0)
  }
  else{
    return(m - indexes[1] + 1)
  }
  # legacy code using a while loop:
  # k <- 0
  # CONT <- TRUE
  # while ((k < m) && CONT) {
  #     if (sorted.pval[k + 1] > lambda/(m - k)) {
  #         CONT <- FALSE
  #     } else {
  #         k <- k + 1
  #     }
  # }
  # return(m - k)
}
# TODO: zeta.HB.sorted that assumes that the pvalues are sorted and doesn't sort them

#' @rdname zeta
#' @export
zeta.trivial <- function(pval, lambda) {
  return(length(pval))
}

#' @rdname zeta
#' @export
zeta.DKWM <- function(pval, lambda) {
  s <- length(pval)
  sorted.pval <- c(0, sort(pval))
  dkwm <- min((sqrt(log(1/lambda)/2)/(2 * (1 - sorted.pval)) + sqrt(log(1/lambda)/(8 * (1 - sorted.pval)^2) + 
                                                                      (s - seq(0, s))/(1 - sorted.pval)))^2,
              na.rm=TRUE)
  return(min(s, floor(dkwm)))
}

# TODO BEFORE MERGE: delete this
nb.elements <- function(C) {
  H <- length(C)
  count <- length(C[[H]])
  if (H > 1) {
    for (h in (H - 1):1) {
      Ch <- C[[h]]
      for (j in 1:length(Ch)) {
        Chj <- Ch[[j]]
        # the following check allows to use an extended
        # tree where leaves are duplicated  at multiple
        # depths so as to appear at the larger depth H
        # or else the function would be simpler : just sum
        # the length of all C[[h]]
        if (Chj[1] < Chj[2])
          count <- count + 1
      }
    }
  }
  return(count)
}

# TODO BEFORE MERGE: rename nb.elements, document
#' Number of unique regions in a reference family with forest structure
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' 
#' @return An integer, the number of regions.
#' 
#' @export
nb.elements.no.extension <- function(C) {
  H <- length(C)
  count <- 0
  for (h in H:1) {
    count <- count + length(C[[h]])
  }
  return(count)
}

# TODO BEFORE MERGE: delete this
zetas.tree <- function(C, leaf_list, method, pvalues, alpha) {
  H <- length(C)
  K <- nb.elements(C)
  leaves <- length(leaf_list)
  zeta_leaves <- numeric(leaves)
  CH <- C[[H]]
  for (i in 1:length(CH)) {
    CHi <- CH[[i]]
    if (CHi[1] == CHi[2]) { # useless leaf check because every bottom region is a leaf? <- NO, not necessarily
      # oh it's worse than that, this piece of code ASSUMES that the tree is extended
      # to have all leaves at the bottom and WILL fail if that's not the case
      # this is a bug
      zeta_leaves[CHi[1]] <- method(pvalues[leaf_list[[CHi[1]]]], alpha/K)
    }
  }
  ZL <- list()
  for (h in H:1) {
    Ch <- C[[h]]
    len <- length(Ch)
    zeta_inter <- numeric(len)
    for (j in 1:len) {
      Chj <- Ch[[j]]
      if (Chj[1] < Chj[2]) {
        pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
        zeta_inter[j] <- method(pvals, alpha/K)
      } else {
        zeta_inter[j] <- zeta_leaves[Chj[1]]
      }
    }
    ZL[[h]] <- zeta_inter
  }
  return(ZL)
}

# TODO BEFORE MERGE: rename zetas.tree, change call of nb.elements
#' Estimate of the proportion of true nulls in each node of a tree
#' 
#' @description
#' Takes a forest structure as input, given by the couple \code{C}, \code{leaf_list} 
#' and returns the corresponding \eqn{\zeta_k}'s according to the method(s) chosen.
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' @param method A function with arguments \code{(pval, lambda)} that can compute an upper bound on the false positives in the region associated to 
#' the \eqn{p}-values \code{pval} at confidence level \code{1 - lambda}. It can also be a list of such functions, where the \code{h}-th function 
#' is used at depth \code{h} in the tree structure, that is on the \eqn{R_k}'s represented by the elements found in \code{C[[h]]}. Finally, it can also be a list 
#' of lists of such functions, mimicking the structure of \code{C} itself, that is, \code{method[[h]][[j]]} is applied the \eqn{R_k} represented by \code{C[[h]][[j]]}.
#' @param pvalues A vector of \eqn{p}-values, must be of size \code{m}, with \code{m} the highest element found in the vectors of \code{leaf_list}.
#' @param alpha A target error level in \eqn{]0,1[]}.
#' @param refine A boolean, \code{FALSE} by default. Whether to use the step-down refinement to try to produce smaller \eqn{\zeta_k}'s, see Details.
#' @param verbose A boolean, \code{FALSE} by default. Whether to print information about the (possibly multiple) round(s) of step-down refinement.
#' 
#' @return \code{ZL}: A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' 
#' @details The proportion of true nulls in each node is estimated by an union bound on the regions. 
#' That is, the provided method(s) is/are applied at level \eqn{\frac{\alpha}{K}} where \eqn{K} is the number of regions. 
#' In the step-down refinement, if we find a \eqn{R_k} with associated \eqn{\zeta_k=0}, that is, we think that the region 
#' contains only false null hypotheses, we can remove it and run again the \eqn{\zeta_k}'s computation using \eqn{K-1} instead of 
#' \eqn{K} in the union bound, and so on until we don't reduce the "effective" number of regions.
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand, G. (2018). Multiple testing and post hoc bounds for heterogeneous data. PhD thesis, see Appendix B.2 for the step-down refinement.
#' @export
#' @examples
#' m <- 1000
#' dd <- dyadic.from.window.size(m, s = 10, method = 2)
#' leaf_list <- dd$leaf_list
#' pvalues <- runif(m)
#' C <- dd$C
#' method <- zeta.trivial
#' ZL <- zetas.tree(C, leaf_list, method, pvalues, alpha = 0.05)
#' ZL
zetas.tree.no.extension <- function(C, leaf_list, method, pvalues, alpha, refine=FALSE, verbose=FALSE) {
  H <- length(C)
  K <- nb.elements.no.extension(C)
  ZL <- list()
  new_K <- K
  continue <- TRUE
  nb_loop <- 0
  while (continue) {
    nb_loop <- nb_loop + 1
    usage_K <- new_K
    new_K <- K
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
        if(typeof(method) == "list"){
          if(typeof(method[[h]]) == "list"){
            zeta_method <- method[[h]][[j]]
          }
          else{
            zeta_method <- method[[h]]
          }
        }
        else{
          zeta_method <- method
        }
        zeta_inter[j] <- zeta_method(pvals, alpha / usage_K)
        if (zeta_inter[j] == 0)
          new_K <- new_K - 1
      }
      ZL[[h]] <- zeta_inter
    }
    if (verbose)
      print(paste0("loop number=", nb_loop, ", usage_K=", usage_K, ", new_K=", new_K))
    continue <- refine && (new_K < usage_K)
  }
  return(ZL)
}

# TODO BEFORE MERGE: delete this
zetas.tree.refined <- function(C, leaf_list, method, pvalues, alpha) {
  H <- length(C)
  K <- nb.elements(C)
  leaves <- length(leaf_list)
  zeta_leaves <- numeric(leaves)
  continue <- TRUE
  new_K <- K
  while (continue) {
    usage_K <- new_K
    new_K <- K
    CH <- C[[H]]
    for (i in 1:length(CH)) {
      CHi <- CH[[i]]
      if (CHi[1] == CHi[2]) {
        pvals <- pvalues[leaf_list[[CHi[1]]]]
        zeta_leaves[CHi[1]] <- method(pvals, alpha/usage_K)
        if (zeta_leaves[CHi[1]] == 0) {
          new_K <- new_K - 1
        }
      }
    }
    ZL <- list()
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        if (Chj[1] < Chj[2]) {
          pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
          zeta_inter[j] <- method(pvals, alpha/usage_K)
          if (zeta_inter[j] == 0) 
            new_K <- new_K - 1
        } else {
          zeta_inter[j] <- zeta_leaves[Chj[1]]
        }
      }
      ZL[[h]] <- zeta_inter
    }
    if (new_K == usage_K) {
      continue <- FALSE
    }
  }
  return(ZL)
}

# TODO BEFORE MERGE: delete this
V.star.all.leaves <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  leaves <- length(leaf_list)
  Vec <- numeric(length = leaves)
  id <- seq_len(length.out = leaves)
  for (i in 1:leaves) {
    # len <- length(intersect(S, leaf_list[[i]]))
    len <- sum(S %in% leaf_list[[i]])
    Vec[i] <- min(ZL[[H]][i], len) # this is the part that assumes that the forest is extended
  }
  if (H > 1) {
    for (h in (H - 1):1) {
      len <- length(C[[h]])
      vec_inter <- numeric(length = leaves)
      new_id <- numeric(len)
      for (j in 1:len) {
        Chj <- C[[h]][[j]]
        # len <- length(intersect(S, leaves))
        lvs <- unlist(leaf_list[Chj[1]:Chj[2]])
        len <- sum(S %in% lvs) # 2 OBJECTS NAMED len !!!!
        rm(lvs)
        ss <- sum(Vec[id[which((Chj[1] <= id) & (Chj[2] >= id))]]) # id is totally useless
        # just use ss <- sum(Vec[Chj[1]:Chj[2]]) lol
        intra <- min(ZL[[h]][j], len, ss)
        vec_inter[Chj[1]] <- intra
        new_id[j] <- Chj[1]
      }
      Vec <- vec_inter
      id <- new_id
    }
  }
  return(sum(Vec[id])) # could also just be sum(Vec) given that Vec has 0 values outside of id
}

# TODO BEFORE MERGE: delete this
V.star.all.leaves.no.id <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves)
  for (i in 1:nb_leaves) {
    len_inter <- sum(S %in% leaf_list[[i]])
    Vec[i] <- min(ZL[[H]][i], len_inter) # this is the part that assumes that the forest is extended
  }
  if (H > 1) {
    for (h in (H - 1):1) {
      nb_regions <- length(C[[h]])
      vec_inter <- numeric(nb_leaves)
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        region_vector <- unlist(leaf_list[Chj[1]:Chj[2]])
        len_inter <- sum(S %in% region_vector)
        sum_succ <- sum(Vec[Chj[1]:Chj[2]]) # this part too assumes that the forest is extended
        res <- min(ZL[[h]][j], len_inter, sum_succ)
        vec_inter[Chj[1]] <- res
      }
      Vec <- vec_inter
    }
  }
  return(sum(Vec))
}

# TODO BEFORE MERGE: rename V.star
#' Post hoc bound on the number of false positives
#' 
#' @description Computes the post hoc upper bound \eqn{V^*(S)} on the number of false positives in a 
#' given selection set \eqn{S} of hypotheses, using a reference family \eqn{(R_k, \zeta_k)} that possess the forest structure
#' (see Reference).
#' 
#' @param S An integer vector with the indices of the hypotheses of the selection set. Does not need to be ordered.
#' @param C A list of list representing the forest structure. Each list of \code{C} represents a level of
#' depth in the forest structure. So, \code{C[[1]]} lists the regions at depth 1, \code{C[[2]]} lists the regions at depth 2,
#' and so on. We then use the fact that each region of the reference family
#' is the union of some of its atom over an integer interval. 
#' So he elements of each list \code{C[[i]]} are integer vectors of length 2 representing this interval. 
#' For example, \code{C[[1]][[1]] = c(1, 5)} means that the first region at depth 1 is the union of the first 5 atoms, 
#' \code{C[[2]][[3]] = c(4, 5)} means that the third region at depth 2 is the union of the atoms 4 and 5, and
#' \code{C[[3]][[5]] = c(5, 5)} means that the fifth region at depth 3 is simply the fifth atom. 
#' @param ZL A list of integer vectors representing the upper bound \eqn{\zeta_k} associated to a region \eqn{R_k} in the reference family. 
#' \code{ZL[[h]][j]} is the \eqn{\zeta_k} associated to the \eqn{R_k} described by \code{C[[h]][[j]]}.
#' @param leaf_list A list of vectors. Each vector is an integer array. The i-th vector contains the indices
#' of the hypotheses in the i-th atom. Atoms form a partition of the set of hypotheses indices : 
#' there cannot be overlap, and each index has to be inside one of the atoms.
#' 
#' @return An integer value, the post hoc upper bound \eqn{V^*(S)}.
#' 
#' @details For \code{V.star}, the forest structure doesn't need to be complete. That is, 
#' in \code{C}, some trivial intervals \code{c(i,i)} corresponding to regions that are atoms may be missing. 
#' 
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @examples
#' m <- 20
#' C <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15), c(17, 17), c(18, 19)),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' ZL <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#' V.star(1, C, ZL, leaf_list)
#' 
#' V.star(1:5, C, ZL, leaf_list)
#' 
#' V.star(13:15, C, ZL, leaf_list)
#' 
#' V.star(1:m, C, ZL, leaf_list)
#' @export
V.star.no.extension <- function(S, C, ZL, leaf_list) {
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
  # this initialization also takes care of the minima
  # between \zeta_k and card(S inter R_k)
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions>0) {
      for (j in 1:nb_regions) {
        Chj <- C[[h]][[j]]
        sum_succ <- sum(Vec[Chj[1]:Chj[2]])
        res <- min(ZL[[h]][j], sum_succ)
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(sum(Vec))
}

#' Prune a forest structure to speed up computations
#' 
#' @description The pruned forest structure makes the computation of [V.star()],
#' [curve.V.star.forest.naive()] and [curve.V.star.forest.fast()] faster.
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' @param prune.leafs A boolean, \code{FALSE} by default. If \code{TRUE}, will also prune atoms/leafs for which \eqn{\zeta_k \geq |R_k|},
#' this makes the computation of [V.star()] and [curve.V.star.forest.naive()] even faster but should not be used with [curve.V.star.forest.fast()]
#' because this needs the structure to be complete (i.e., with all its atoms). This is why the default option is \code{FALSE}.
#' 
#' @return A list with three named elements. \describe{
#' \item{\code{VstarNm}}{\eqn{V^*(\mathbb N_m)} is computed as by-product by the algorithm, so we might as well return it.}
#' \item{\code{C}}{The new \code{C} after pruning.}
#' \item{\code{ZL}}{The new \code{ZL} after pruning.}
#' }
#' 
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G., preprint to appear with the description of pruning
#' @export
pruning <- function(C, ZL, leaf_list, prune.leafs = FALSE) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  Vec <- numeric(nb_leaves) 
  for (i in 1:nb_leaves) {
    Vec[i] <- length(leaf_list[[i]])
    # The initialization term for each atom P_i
    # is equivalent to completing the family if it isn't,
    # assuming that leaf_list does indeed contain all leaves
    # and some were just eventually missing in C and ZL.
    # Using the length of the atoms assure that we 
    # will catch atoms with a \zeta_i larger than its length
  }
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions > 0) {
      for (j in nb_regions:1) {
        Chj <- C[[h]][[j]]
        candidate <- ZL[[h]][j]
        sum_succ <- sum(Vec[Chj[1]:Chj[2]])
        if (candidate >= sum_succ) {
          res <- sum_succ
          if(prune.leafs || (Chj[1] < Chj[2])){
            # Either the region is not a leaf, 
            # or it's a leaf and we chose to prune the leafs
            # so we prune.
            ZL[[h]] <- ZL[[h]][-j]
            C[[h]][[j]] <- NULL
          }
          else{
            # The region is a leaf/atom and we chose
            # to keep the leafs, so we keep it,
            # but we change the \zeta_i to the length of P_i.
            ZL[[h]][j] <- sum_succ
          }
        } else {
          res <- candidate
        }
        Vec[Chj[1]:Chj[2]] <- 0
        Vec[Chj[1]] <- res
      }
    }
  }
  return(list(VstarNm = sum(Vec),
              C = C,
              ZL = ZL
  )
  )
}
# TODO later: fast curve.Vmax is slower with the pruning
# which is unexpected, maybe it is because the pruning leaves gaps
# => maybe revamp to delete gaps? 
# (??? curve.Vmax not found, maybe this TODO should be deleted)

#' Compute a curve of post hoc bounds based on a reference family with forest structure
#' 
#' @name curve.V.star.forest
#' 
#' @description
#' Computes the post hoc upper bound \eqn{V^*(S_t)} on the number of false positives in a 
#' given sequence of selection sets \eqn{S_t} of hypotheses, such that
#' \eqn{S_t\subset S_{t+1}} and \eqn{|S_t| = t}, 
#' using a reference family \eqn{(R_k, \zeta_k)} that possess the forest structure
#' (see References).
#' 
#' @details Two functions are available
#' \describe{
#' \item{\code{curve.V.star.forest.naive}}{Repeatedly calls [V.star()] on each \eqn{S_t}, which is not optimized and time-consuming, this should be used in practice.}
#' \item{\code{curve.V.star.forest.fast}}{A fast and optimized version that leverage the fact that\eqn{S_{t+1}} is the union of \eqn{S_t} and a single hypothesis index. This 
#' option first completes the forest, because the algorithm needs that, and the completion fails if the input is a pruned forest (see [pruning()]), 
#' so if a pruned forest is given as input, it MUST be said with the \code{is.pruned} argument 
#' so that the function skips completion.}
#' }
#'
#' @param perm An integer vector of elements in \code{1:m}, all different, and of size up to \code{m} (in which case it's a permutation, hence the name). 
#' The set \eqn{S_t} is represented by \code{perm[1:t]}.
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' @param pruning A boolean, \code{FALSE} by default. Whether to prune the forest (see [pruning()]) before computing the bounds.
#' @param is.pruned A boolean, \code{FALSE} by default. If \code{TRUE}, assumes that the forest structure has already been completed and then pruned 
#' and so skips the completion step. Must be set to \code{TRUE} if giving a pruned forest because in that case the completion step must be skipped or it will fail.
#' 
#' @return A vector of length of same length as \code{perm}, where the \code{t}-th
#' element is \eqn{V^*(S_t)}.
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G., preprint to appear with the description of pruning and of the fast algorithm to compute the curve.
#' @export
#' @examples
#' m <- 20
#' C <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(c(4, 5), c(8, 9), c(10, 10), c(12, 12), c(13, 15), c(17, 17), c(18, 19)),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' ZL <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#' curve.V.star.forest.naive(1:m, C, ZL, leaf_list, pruning = FALSE)
#' 
#' curve.V.star.forest.naive(1:m, C, ZL, leaf_list, pruning = TRUE)
#' 
#' curve.V.star.forest.fast(1:m, C, ZL, leaf_list, pruning = FALSE)
#' 
#' curve.V.star.forest.fast(1:m, C, ZL, leaf_list, pruning = TRUE)


# TODO BEFORE MERGE: change call of V.star
#' @rdname curve.V.star.forest
#' @export
curve.V.star.forest.naive <- function(perm, C, ZL, leaf_list, pruning = FALSE){
  vstars <- numeric(length(perm))
  
  # the naive version doesn't need a proper completion of the
  # forest structure because V.star.no.extension
  # implicitly completes, and for the same reason
  # it can use super pruning
  
  if (pruning){
    pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE)
    C <- pruned$C
    ZL <- pruned$ZL
    m <- length(unlist(leaf_list))
    if(length(perm) == m){
      # means that the last bound to compute
      # is V^*({1, ..., m}),
      # but the pruning already computed
      # V^*({1, ..., m}) as a by-product so we
      # might as well use it:
      vstars[m] <- pruned$VstarNm
      perm <- perm[-m]
    }
  }
  
  S <- numeric(0)
  j <- 0
  for (t in perm){
    j <- j + 1
    S <- c(S, t)
    vstars[j] <- V.star.no.extension(S, C, ZL, leaf_list) 
  }
  return(vstars)
}

#' @rdname curve.V.star.forest
#' @export
curve.V.star.forest.fast <- function(perm, C, ZL, leaf_list, pruning = FALSE, is.pruned = FALSE){
  
  vstars <- numeric(length(perm))
  
  if (! is.pruned) {
    # the fast version needs a proper completion of the
    # forest structure, and for the same reason
    # it must not use super pruning
    completed <- forest.completion(C, ZL, leaf_list)
    C <- completed$C
    ZL <- completed$ZL
    
    if (pruning) {
      is.pruned <- TRUE # useless atm, but kept for clarity
      pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE)
      C <- pruned$C
      ZL <- pruned$ZL
      m <- length(unlist(leaf_list))
      
      if (length(perm) == m) {
        # means that the last bound to compute
        # is V^*({1, ..., m}),
        # but the pruning already computed
        # V^*({1, ..., m}) as a by-product so we
        # might as well use it:
        vstars[m] <- pruned$VstarNm
        perm <- perm[-m]
      }
    }
  }
  
  H <- length(C)
  
  etas <- ZL
  K.minus <- list()
  for (h in 1:H){
    etas[[h]] <- rep(0, length(ZL[[h]]))
    K.minus[[h]] <- list()
    if (length(ZL[[h]]) > 0){
      for (j in 1:length(ZL[[h]])){
        if (ZL[[h]][j] == 0){
          K.minus[[h]][[j]] <- C[[h]][[j]]
        }
      }
    }
  }
  
  for (t in 1:length(perm)) {
    
    i.t <- perm[t]
    if (t > 1) {
      previous.vstar <- vstars[t - 1]
    } else {
      previous.vstar <- 0
    }
    
    ################################
    # SEARCHING IF i_t IS IN K MINUS
    # if so, go.next == TRUE
    # and we just go next to step t+1
    go.next <- FALSE
    for (h in 1:H) {
      if (go.next) {
        break
      }
      for (couple in K.minus[[h]]) {
        if (! is.null(couple)) {
          lower_leaf <- leaf_list[[couple[1]]]
          lower_hyp <- lower_leaf[1]
          upper_leaf <- leaf_list[[couple[2]]]
          upper_hyp <- upper_leaf[length(upper_leaf)]
          if ((i.t >= lower_hyp) && (i.t <= upper_hyp)) {
            go.next <- TRUE
            # print(paste0(i.t, " is in K minus"))
            break
          }
        }
      }
    }
    # print(paste0(i.t, " isn't in K minus"))
    #########################################
    
    # COMPUTING V.STAR AND UPDATING K.MINUS AND ETAS
    ################################################
    if (go.next) {
      vstars[t] <- previous.vstar
    } else {
      # Here, i_t isn't in K minus
      for (h in 1:H) {
        nb_regions <- length(C[[h]])
        if(nb_regions > 0){
          is.found <- FALSE
          for (j in 1:nb_regions) {
            couple <- C[[h]][[j]]
            lower_leaf <- leaf_list[[couple[1]]]
            lower_hyp <- lower_leaf[1]
            upper_leaf <- leaf_list[[couple[2]]]
            upper_hyp <- upper_leaf[length(upper_leaf)]
            if((i.t >= lower_hyp) && (i.t <= upper_hyp)){
              # we found k^{(t,h)}
              is.found <- TRUE
              break
            }
          }
          if (! is.found) {
            next
          }
          etas[[h]][[j]] <- etas[[h]][[j]] + 1
          if(etas[[h]][[j]] < ZL[[h]][[j]]){
            # pass
          } else {
            K.minus[[h]][[j]] <- C[[h]][[j]]
            break
          }
        }
      }
      vstars[t] <- previous.vstar + 1
    }
    ################################################
    
  }
  return(vstars)
}

# the forest must not be pruned beforehand
# the following code fails if the input is a pruned forest
# TODO BEFORE MERGE: document
#' Complete a forest structure
#' 
#' @description Completes the forest in the sens of the Reference: adds the missing atoms/leafs 
#' in the reference family with a forest structure \eqn{(R_k, \zeta_k)} 
#' so that each atom is well represented by a \eqn{R_k}. The associated \eqn{\zeta_k} is 
#' taken as the trivial \eqn{|R_k|}.
#' 
#' @details The forest must not be pruned (with [pruning()]) beforehand. The code will not behave expectedly 
#' and will return a wrong result if a pruned forest is given as input. Maybe the function could be rewritten 
#' going from the leaves to the roots instead of the contrary, to avoid this issue.
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' 
#' @return A list with two named elements. \describe{
#' \item{\code{C}}{The new \code{C} after completion.}
#' \item{\code{ZL}}{The new \code{ZL} after completion.}
#' }
#' 
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @export
forest.completion <- function(C, ZL, leaf_list) {
  H <- length(C)
  
  leaves.to.place <- 1:length(leaf_list)
  len.to.place <- length(leaves.to.place)
  to.delete <- numeric(0)
  
  for (h in 1:H) {
    
    j <- 1
    l <- 1
    
    while (j <= len.to.place) {
      expected_leaf <- leaves.to.place[j]
      end_of_line <- l > length(C[[h]])
      if (! end_of_line) {Chl <- C[[h]][[l]]}
      if ((! end_of_line) && expected_leaf == Chl[[1]]) {
        if (expected_leaf == Chl[[2]]) {
          to.delete <- c(to.delete, j)
        }
        j <- j + Chl[2] - Chl[1] +1
      } else {
        C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)), l - 1)
        ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]), l - 1)
        to.delete <- c(to.delete, j)
        j <- j + 1
      }
      l <- l + 1
    }
    
    leaves.to.place <- leaves.to.place[-to.delete]
    len.to.place <- length(leaves.to.place)
    to.delete <- numeric(0)
    
  }
  
  if (len.to.place > 0) {
    h <- H + 1
    C[[h]] <- list()
    ZL[[h]] <- numeric(0)
    for (expected_leaf in leaves.to.place) {
      C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)))
      ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]))
    }
  }
  
  return(list(C = C, ZL = ZL))
}

# TODO BEFORE MERGE: delete this
V.star <- function(S, C, ZL, leaf_list) {
  all_leaves <- tree.expand(C, ZL, leaf_list)
  return(V.star.all.leaves(S, 
                           all_leaves$C, 
                           all_leaves$ZL, 
                           leaf_list))
}

# TODO BEFORE MERGE: delete this
V.star.no.id <- function(S, C, ZL, leaf_list) {
  all_leaves <- tree.expand(C, ZL, leaf_list)
  return(V.star.all.leaves.no.id(S, 
                                 all_leaves$C, 
                                 all_leaves$ZL, 
                                 leaf_list))
}

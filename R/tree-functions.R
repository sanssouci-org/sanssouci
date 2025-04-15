#' Create a complete dyadic tree structure
#' 
#' @description
#' Produce a set of regions with forest structure, with a single tree, which is dyadic.
#' 
#' 
#' @name dyadic
#' @details \describe{
#' \item{\code{dyadic.from.leaf_list}}{Dyadic tree structure from a given list of atoms/leafs}
#' \item{\code{dyadic.from.window.size}}{Dyadic tree structure from window size: the number of elements in each atom/leaf is set to \code{s}}
#' \item{\code{dyadic.from.height}}{Dyadic tree structure from height: the total height of the tree is set to \code{H}}
#' }
#' 
#' @param method A numeric value. If \code{method == 1}, start from the leaves
#' and group nodes of a same height 2 by 2 as long as possible. If
#' \code{method==2}, start from the root and divide nodes in 2 nodes of equal
#' size as long as possible
#' 
#' @return A list with two named elements:\describe{
#' \item{\code{leaf_list}}{A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.}
#' \item{\code{C}}{A list of list representing the forest structure. See [V.star()] for more information.}
#' }
#' 
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
#' @examples
#' m <- 6
#' dd <- dyadic.from.window.size(m, s = 2, method = 2)
#' str(dd)
#' 
#' dd <- dyadic.from.height(m, H = 3, method = 2)
#' str(dd)
#' 
#' dd <- dyadic.from.height(m, method = 2)
#' str(dd)
#'
#' leaf_list <- dd$leaf_list
#' dd <- dyadic.from.leaf_list(leaf_list, method = 2)
#' str(dd)
NULL

#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
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
  return(list(leaf_list = leaf_list, C = C))
}

#' @param m An integer value, the number of hypotheses to have in the structure
#' @param s An integer value, the number of elements in each atom/leaf
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
  C <- dyadic.from.leaf_list(leaf_list, method)$C
  return(list(leaf_list = leaf_list, C = C))
}

#' @param H An integer value, the desired maximal height of the tree. If NULL (by default), 
#' use \code{floor(2 + log2(m - 1))} which gives the maximum achievable height given \code{m}.
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
  C <- dyadic.from.leaf_list(leaf_list, method)$C
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
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
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
  if (! length(indexes)){
    return(0)
  }
  else{
    return(m - indexes[1] + 1)
  }
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
  dkwm <- min((sqrt(log(1/lambda)/2)/(2 * (1 - sorted.pval)) + 
                  sqrt(log(1/lambda)/(8 * (1 - sorted.pval)^2) + 
                         (s - seq(0, s))/(1 - sorted.pval)
                       )
                )^2,
              na.rm=TRUE)
  return(min(s, floor(dkwm)))
}

#' Number of unique regions in a reference family with forest structure
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' 
#' @return An integer, the number of regions.
#' 
#' @export
nb.elements <- function(C) {
  H <- length(C)
  count <- 0
  for (h in H:1) {
    count <- count + length(C[[h]])
  }
  return(count)
}

# TODO : allow zetas.tree to take additional arguments to pass to method

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
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
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
zetas.tree <- function(C, leaf_list, method, pvalues, alpha, refine=FALSE, verbose=FALSE) {
  H <- length(C)
  K <- nb.elements(C)
  ZL <- list()
  new_K <- K
  continue <- TRUE
  nb_loop <- 0
  while (continue) {
    usage_K <- new_K
    new_K <- K
    for (h in H:1) {
      Ch <- C[[h]]
      len <- length(Ch)
      zeta_inter <- numeric(len)
      for (j in 1:len) {
        Chj <- Ch[[j]]
        pvals <- pvalues[unlist(leaf_list[Chj[1]:Chj[2]])]
        if(typeof(method) == "list") {
          if(typeof(method[[h]]) == "list") {
            zeta_method <- method[[h]][[j]]
          }
          else {
            zeta_method <- method[[h]]
          }
        }
        else {
          zeta_method <- method
        }
        zeta_inter[j] <- zeta_method(pvals, alpha / usage_K)
        if (zeta_inter[j] == 0)
          new_K <- new_K - 1
      }
      ZL[[h]] <- zeta_inter
    }
    if (verbose) {
      nb_loop <- nb_loop + 1
      print(paste0("loop number=", nb_loop, ", usage_K=", usage_K, ", new_K=", new_K))
    }
    continue <- refine && (new_K < usage_K)
  }
  return(ZL)
}

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
V.star <- function(S, C, ZL, leaf_list) {
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
#' @param delete.gaps A boolean, \code{FALSE} by default. If \code{TRUE}, will also delete the gaps in the structure induced 
#' by the pruning, see [delete.gaps()].
#' 
#' @return A list with three named elements. \describe{
#' \item{\code{VstarNm}}{\eqn{V^*(\mathbb N_m)} is computed as by-product by the algorithm, so we might as well return it.}
#' \item{\code{C}}{The new \code{C} after pruning.}
#' \item{\code{ZL}}{The new \code{ZL} after pruning.}
#' }
#' 
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
#' @export
pruning <- function(C, ZL, leaf_list, prune.leafs = FALSE, delete.gaps = FALSE) {
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
  if (delete.gaps) {
    gaps.deleted <- delete.gaps(C, ZL, leaf_list)
    C <- gaps.deleted$C
    ZL <- gaps.deleted$ZL
  }
  return(list(VstarNm = sum(Vec),
              C = C,
              ZL = ZL
  )
  )
}

#' Delete the gaps induced by pruning
#' 
#' @description
#' A small optimization that can be done after pruning
#' that can speed up computations (it removes the gaps introduced by the pruning)
#' 
#' @details
#' See [pruning()].
#' 
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' 
#' @return A list with two named elements. \describe{
#' \item{\code{C}}{The new \code{C} after deleting the gaps.}
#' \item{\code{ZL}}{The new \code{ZL} after deleting the gaps.}
#' }
#' 
#' @export
delete.gaps <- function(C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  continue <- TRUE
  newC <- list()
  newZL <- list()
  loop.counter <- 1
  while(continue) {
    newC[[loop.counter]] <- list()
    newZL[[loop.counter]] <- numeric(0)
    leaf.available <- ! logical(nb_leaves)
    regions.delete <- list()
    for (h in 1:H) {
      nb_regions <- length(C[[h]])
      if (nb_regions > 0) {
        for (l in 1:nb_regions) {
          Chl <- C[[h]][[l]]
          if (all(leaf.available[Chl[1]:Chl[2]])) {
            newC[[loop.counter]][[Chl[1]]] <- Chl
            newZL[[loop.counter]][Chl[1]]  <- ZL[[h]][[l]]
            leaf.available[Chl[1]:Chl[2]] <- FALSE
            regions.delete <- append(regions.delete, list(c(h, l)))
          }
        }
      }
    }
    for (couple in rev(regions.delete)) {
      h <- couple[1]
      l <- couple[2]
      C[[h]][[l]] <- NULL
      ZL[[h]] <- ZL[[h]][-l]
    }
    if(length(newC[[loop.counter]]) > 0) {
      for(l in length(newC[[loop.counter]]):1) {
        if (is.null(newC[[loop.counter]][[l]])) {
          newC[[loop.counter]][[l]] <- NULL
          newZL[[loop.counter]] <- newZL[[loop.counter]][-l]
        }
      }
    }
    loop.counter <- loop.counter + 1
    continue <- any(sapply(X = ZL, FUN = function(vec){length(vec) > 0}))
  }
  return(list(C = newC, ZL = newZL))
}

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
#' \item{\code{curve.V.star.forest.fast}}{A fast and optimized version that leverage the fact that\eqn{S_{t+1}} is the union of \eqn{S_t} and a single hypothesis index. 
#' The algorithm needs to work on a complete forest, so this version first completes the forest (unless told that the forest has already been completed, see [forest.completion()]), 
#' and the completion fails if the input is a pruned forest (see [pruning()]), 
#' so if a pruned forest is given as input, it MUST be said with the \code{is.pruned} argument 
#' so that the function skips completion (so the pruned forest given as input must also be complete).}
#' }
#'
#' @param perm An integer vector of elements in \code{1:m}, all different, and of size up to \code{m} (in which case it's a permutation, hence the name). 
#' The set \eqn{S_t} is represented by \code{perm[1:t]}.
#' @param C A list of list representing the forest structure. See [V.star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds \eqn{\zeta_k} of the forest structure. See [V.star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest structure. See [V.star()] for more information.
#' @param pruning A boolean, \code{FALSE} by default. Whether to prune the forest (see [pruning()]) before computing the bounds. Ignored if \code{is.pruned} is \code{TRUE}.
#' @param is.pruned A boolean, \code{FALSE} by default. If \code{TRUE}, assumes that the forest structure has already been completed (see [forest.completion()]) and then pruned (see [pruning()])
#' and so skips the completion step and optional pruning step. Must be set to \code{TRUE} if giving a pruned forest, see Details.
#' @param is.complete A boolean, \code{FALSE} by default. If \code{TRUE}, assumes that the forest structure has already been completed (see [forest.completion()]) and so skips the completion step.
#' Ignored if \code{is.pruned} is \code{TRUE}.
#' @param delete.gaps A boolean, \code{FALSE} by default. If \code{TRUE}, will also delete the gaps in the structure induced 
#' by the pruning, see [delete.gaps()]. Ignored if \code{pruning} is \code{FALSE}.
#' 
#' @return A vector of length of same length as \code{perm}, where the \code{t}-th
#' element is \eqn{V^*(S_t)}.
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
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

#' @rdname curve.V.star.forest
#' @export
curve.V.star.forest.naive <- function(perm, C, ZL, leaf_list, pruning = FALSE, delete.gaps = FALSE){
  vstars <- numeric(length(perm))
  
  # the naive version doesn't need a proper completion of the
  # forest structure because V.star.no.extension
  # implicitly completes, and for the same reason
  # it can use super pruning
  
  if (pruning){
    pruned <- pruning(C, ZL, leaf_list, prune.leafs = TRUE, delete.gaps = delete.gaps)
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
    vstars[j] <- V.star(S, C, ZL, leaf_list) 
  }
  return(vstars)
}

#' @rdname curve.V.star.forest
#' @export
curve.V.star.forest.fast <- function(perm, C, ZL, leaf_list, pruning = FALSE, is.pruned = FALSE, is.complete = FALSE, delete.gaps = FALSE){
  
  vstars <- numeric(length(perm))
  
  if (! is.pruned) {
    
    if (! is.complete) {
      # the fast version needs a proper completion of the
      # forest structure, and for the same reason
      # it must not prune the leaves
      completed <- forest.completion(C, ZL, leaf_list)
      C <- completed$C
      ZL <- completed$ZL
    }
    
    if (pruning) {
      is.pruned <- TRUE # useless atm, but kept for clarity
      pruned <- pruning(C, ZL, leaf_list, prune.leafs = FALSE, delete.gaps = delete.gaps)
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
  m <- length(unlist(leaf_list))
  
  # preparation of the initial value of the etas (a copy of zetas but 
  # with only zeroes), of the initial value of K^- (with the R_k's
  # such that zeta_k = 0) and of a m x H matrix M such that
  # M[i, h] gives the index k of C[[h]] such that hypothesis i
  # is in the R_k given by C[[h]][[k]]
  etas <- ZL
  K.minus <- list()
  M <- matrix(0, ncol = H, nrow = m)
  for (h in 1:H){
    zeta_depth_h <- ZL[[h]]
    length_zeta_depth_h <- length(zeta_depth_h)
    etas[[h]] <- rep(0, length_zeta_depth_h)
    K.minus[[h]] <- vector("list", length(C[[h]]))
    if (length_zeta_depth_h > 0){
      for (k in 1:length_zeta_depth_h){
        if (zeta_depth_h[k] == 0){
          K.minus[[h]][[k]] <- C[[h]][[k]]
        }
        first_leaf <- leaf_list[[C[[h]][[k]][1]]]
        last_leaf <- leaf_list[[C[[h]][[k]][2]]]
        M[first_leaf[1]:last_leaf[length(last_leaf)], h] <- k
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
    
    # SEARCHING IF i_t IS IN K MINUS
    # if so, go.next == TRUE
    # and we just go next to step t+1
    go.next <- FALSE
    for (h in 1:H) {
      k <- M[i.t, h]
      if ((k > 0) && (! is.null(K.minus[[h]][[k]]))){
        go.next <- TRUE
        break
      }
    }
    # print(paste0(i.t, " isn't in K minus"))
    
    # COMPUTING V.STAR AND UPDATING K.MINUS AND ETAS
    if (go.next) {
      vstars[t] <- previous.vstar
    } else {
      # Here, i_t isn't in K minus
      for (h in 1:H) {
        k <- M[i.t, h]
        if (k > 0){
          # if k == 0,
          # there is no k^{(t,h)} because there is a 
          # gap in the structure (because of pruning)
          # in this case we don't do anything
          etas[[h]][[k]] <- etas[[h]][[k]] + 1
          if(etas[[h]][[k]] >= ZL[[h]][[k]]){
            K.minus[[h]][[k]] <- C[[h]][[k]]
            break
          }
        }
      }
      vstars[t] <- previous.vstar + 1
    }
    
  }
  return(vstars)
}

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
#' @references Durand G. (2025). A fast algorithm to compute a curve of confidence upper bounds for the False Discovery Proportion using a reference family with a forest structure. arXiv:2502.03849.
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

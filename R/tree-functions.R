#' Estimate the number of true null hypotheses among a set of p-values
#'
#' @description
#' An upper bound on the number of true null hypotheses in the region associated
#' to the \eqn{p}-values \code{p_values}
#' is computed with confidence \code{1 - lambda}.
#' The functions described here can be used as the \code{method} argument
#' of [zetas_tree()].
#'
#' @param p_values A vector of \eqn{p}-values.
#' @param lambda A numeric value in \eqn{[0,1]},
#' the target error level of the test.
#' @param k The positive integer \eqn{k} used by the
#' \eqn{k}-Bonferroni procedure. By default it is equal to 1.
#' @param ... Additional arguments that may be passed
#' to specific \code{zeta} functions.
#' @name zeta
#' @return The number of true nulls is over-estimated as follows:
#' \describe{
#' \item{\code{zeta_dkwm}}{
#' Inversion of the Dvoretzky-Kiefer-Wolfowitz-Massart inequality
#' (related to the Storey estimator of the proportion of true nulls)
#' with parameter \code{lambda}
#' }
#' \item{\code{zeta_hb}}{
#' Number of conserved hypotheses by the Holm-Bonferroni
#' procedure with parameter \code{lambda}
#' }
#' \item{\code{zeta_kbonf}}{
#' Number of conserved hypotheses by the \eqn{k}-Bonferroni
#' procedure with parameter \code{lambda}, plus \eqn{k-1}
#' }
#' \item{\code{zeta_trivial}}{
#' The size of the p-value set which is the trivial
#' upper bound (\code{lambda} is not used)
#' }
#' }
#' @details The \eqn{k}-Bonferroni procedure controls the
#' \eqn{k}-Familywise Error Rate (FWER) at the desired level,
#' hence the number of conserved hypotheses, plus \eqn{k-1},
#' is a suitable upper bound
#' (because up to \eqn{k-1} rejected hypotheses are also true nulls).
#' For \eqn{k=1} (the default), it is the regular Bonferroni procedure.
#' The Holm-Bonferroni procedure controls the FWER, hence the number
#' of conserved hypotheses is a suitable upper bound.
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956).
#' Asymptotic minimax character of the sample distribution function and of
#' the classical multinomial estimator.
#' The Annals of Mathematical Statistics, pages 642-669.
#' @references Holm, S. A simple sequentially rejective multiple test procedure.
#' Scandinavian Journal of Statistics 6 (1979), pp. 65-70.
#' @references Massart, P. (1990). The tight constant in the
#' Dvoretzky-Kiefer-Wolfowitz inequality.
#' The Annals of Probability, pages 1269-1283.
#' @references Storey, J. D. (2002). A direct approach to
#' false discovery rates.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#' 64(3):479-498.
#' @references Durand G. (2025). Fast confidence bounds for the false
#' discovery proportion over a path of hypotheses. Computo.
#' @examples
#' x <- rnorm(100, mean = c(rep(c(0, 2), each = 50)))
#' p_values <- 1 - pnorm(x)
#' lambda <- 0.05
#' zeta_trivial(p_values, lambda)
#'
#' zeta_hb(p_values, lambda)
#'
#' zeta_dkwm(p_values, lambda)
NULL

#' @rdname zeta
#' @export
zeta_hb <- function(p_values, lambda, ...) {
  m <- length(p_values)
  sorted_p_values <- sort(p_values)

  thresholds <- lambda / (m - 1:m + 1)
  v <- sorted_p_values - thresholds
  indexes <- which(v > 0)
  if (!length(indexes)) {
    return(0)
  } else {
    return(m - indexes[1] + 1)
  }
}
# TODO: zeta_hb.sorted that assumes that the p_values
# are sorted and doesn't sort them

#' @rdname zeta
#' @export
zeta_kbonf <- function(p_values, lambda, k = 1, ...) {
  m <- length(p_values)
  threshold <- lambda * k / m
  indexes <- which(p_values <= threshold)
  bound <- m - length(indexes) + (k - 1)
  return(min(bound, m))
}

#' @rdname zeta
#' @export
zeta_trivial <- function(p_values, lambda, ...) {
  return(length(p_values))
}

#' @rdname zeta
#' @export
zeta_dkwm <- function(p_values, lambda, ...) {
  s <- length(p_values)
  sorted_p_values <- c(0, sort(p_values))
  dkwm <- min(
    (
      sqrt(log(1 / lambda) / 2) / (2 * (1 - sorted_p_values)) +
        sqrt(log(1 / lambda) / (8 * (1 - sorted_p_values)^2) +
               (s - seq(0, s)) / (1 - sorted_p_values))
    )^2,
    na.rm = TRUE
  )
  return(min(s, floor(dkwm)))
}

#' Number of unique regions in a reference family with forest structure
#'
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#'
#' @return An integer, the number of regions.
#'
#' @export
nb_elements <- function(forest_structure) {
  H <- length(forest_structure)
  count <- 0
  for (h in H:1) {
    count <- count + length(forest_structure[[h]])
  }
  return(count)
}

#' Estimate of the proportion of true nulls in each node of a tree
#'
#' @description
#' Takes a forest structure as input, given by
#' the couple \code{C}, \code{leaf_list}
#' and returns the corresponding \eqn{\zeta_k}'s
#' according to the method(s) chosen.
#'
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest
#' structure. See [v_star()] for more information.
#' @param method A function with arguments \code{(p_values, lambda, ...)}
#' that can
#' compute an upper bound on the false positives in the region associated to
#' the \eqn{p}-values \code{p_values} at confidence level \code{1 - lambda}.
#' It can also be a list of such functions, where the \code{h}-th function
#' is used at depth \code{h} in the tree structure, that is on the \eqn{R_k}'s
#' represented by the elements found in \code{C[[h]]}.
#' Finally, it can also be a list
#' of lists of such functions, mimicking the structure of \code{C} itself,
#' that is, \code{method[[h]][[j]]} is applied the \eqn{R_k}
#' represented by \code{C[[h]][[j]]}.
#' @param p_values A vector of \eqn{p}-values, must be of size \code{m},
#' with \code{m} the highest element found in the vectors of \code{leaf_list}.
#' @param alpha A target error level in \eqn{]0,1[]}.
#' @param refine A boolean, \code{FALSE} by default.
#' Whether to use the step-down refinement to try to
#' produce smaller \eqn{\zeta_k}'s, see Details.
#' @param verbose A boolean, \code{FALSE} by default.
#' Whether to print information about the (possibly multiple)
#' round(s) of step-down refinement.
#' @param ... Additional arguments that may be passed to
#' specific \code{zeta} functions.
#'
#' @return \code{ZL}: A list of integer vectors representing
#' the upper bounds \eqn{\zeta_k} of the forest structure.
#' See [v_star()] for more information.
#'
#' @details The proportion of true nulls in each node is estimated by
#' an union bound on the regions.
#' That is, the provided method(s) is/are applied at level
#' \eqn{\frac{\alpha}{K}} where \eqn{K} is the number of regions.
#' In the step-down refinement, if we find a \eqn{R_k} with associated
#' \eqn{\zeta_k=0}, that is, we think that the region
#' contains only false null hypotheses, we can remove it and run again
#' the \eqn{\zeta_k}'s computation using \eqn{K-1} instead of
#' \eqn{K} in the union bound, and so on until we don't reduce
#' the "effective" number of regions.
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand, G. (2018). Multiple testing and
#' post hoc bounds for heterogeneous data.
#' PhD thesis, see Appendix B.2 for the step-down refinement.
#' @references Durand G. (2025). Fast confidence bounds for the false discovery
#' proportion over a path of hypotheses. Computo.
#' @export
#' @examples
#' m <- 1000
#' dd <- dyadic_from_window_size(m, s = 10, method = 2)
#' leaf_list <- dd$leaf_list
#' p_values <- runif(m)
#' C <- dd$C
#' method <- zeta_trivial
#' ZL <- zetas_tree(C, leaf_list, method, p_values, alpha = 0.05)
#' ZL
zetas_tree <- function(C,
                       leaf_list,
                       method,
                       p_values,
                       alpha,
                       refine = FALSE,
                       verbose = FALSE,
                       ...) {
  H <- length(C)
  K <- nb_elements(C)
  ZL <- list()
  NEW_K <- K
  continue <- TRUE
  nb_loop <- 0
  while (continue) {
    USAGE_K <- NEW_K
    NEW_K <- K
    for (h in H:1) {
      C_H <- C[[h]]
      len <- length(C_H)
      zeta_inter <- numeric(len)
      if (len > 0) {
        for (k in 1:len) {
          C_HK <- C_H[[k]]
          extracted_p_values <- p_values[unlist(leaf_list[C_HK[1]:C_HK[2]])]
          args_zeta <- list(...)
          if ("p_cdfs" %in% names(args_zeta)) {
            args_zeta$p_cdfs <- args_zeta$p_cdfs[unlist(
              leaf_list[C_HK[1]:C_HK[2]]
            )]
          }
          if (typeof(method) == "list") {
            if (typeof(method[[h]]) == "list") {
              zeta_method <- method[[h]][[k]]
            } else {
              zeta_method <- method[[h]]
            }
          } else {
            zeta_method <- method
          }
          zeta_inter[k] <- do.call(
            zeta_method,
            c(
              list(
                p_values = extracted_p_values,
                lambda = alpha / USAGE_K
              ),
              args_zeta
            )
          )
          if (zeta_inter[k] == 0) {
            NEW_K <- NEW_K - 1
          }
        }
      }
      ZL[[h]] <- zeta_inter
    }
    if (verbose) {
      nb_loop <- nb_loop + 1
      print(paste0(
        "loop number=",
        nb_loop,
        ", USAGE_K=",
        USAGE_K,
        ", NEW_K=",
        NEW_K
      ))
    }
    continue <- refine && (NEW_K < USAGE_K)
  }
  return(ZL)
}

#' Post hoc bound on the number of false positives
#'
#' @description Computes the post hoc upper bound \eqn{V^*(S)} on the number
#' of false positives in a
#' given selection set \eqn{S} of hypotheses,
#' using a reference family \eqn{(R_k, \zeta_k)}
#' that possess the forest structure
#' (see Reference).
#'
#' @param S An integer vector with the indices of the hypotheses of
#' the selection set. Does not need to be ordered.
#' @param C A list of list representing the forest structure.
#' Each list of \code{C} represents a level of
#' depth in the forest structure.
#' So, \code{C[[1]]} lists the regions at depth 1,
#' \code{C[[2]]} lists the regions at depth 2,
#' and so on. We then use the fact that each region of the reference family
#' is the union of some of its atom over an integer interval.
#' So he elements of each list \code{C[[i]]} are integer vectors of length 2
#' representing this interval.
#' For example, \code{C[[1]][[1]] = c(1, 5)} means that the first region
#' at depth 1 is the union of the first 5 atoms,
#' \code{C[[2]][[3]] = c(4, 5)} means that the third region at depth 2
#' is the union of the atoms 4 and 5, and
#' \code{C[[3]][[5]] = c(5, 5)} means that the fifth region at depth 3
#' is simply the fifth atom.
#' @param ZL A list of integer vectors representing the upper bound
#' \eqn{\zeta_k} associated to a region \eqn{R_k} in the reference family.
#' \code{ZL[[h]][j]} is the \eqn{\zeta_k} associated to the
#' \eqn{R_k} described by \code{C[[h]][[j]]}.
#' @param leaf_list A list of vectors. Each vector is an integer array.
#' The i-th vector contains the indices
#' of the hypotheses in the i-th atom.
#' Atoms form a partition of the set of hypotheses indices :
#' there cannot be overlap, and each index has to be inside one of the atoms.
#' It is really important, for all functions
#' using \code{leaf_list} in this package, that each vector is sorted.
#' Because functions of the package
#' implicitly use that
#' \code{leaf_list[[i]][1]} is the lowest hypothesis
#' index in leaf \eqn{P_i} and that
#' \code{leaf_list[[i]][length(leaf_list[[i]])]} is the largest
#'
#' @return An integer value, the post hoc upper bound \eqn{V^*(S)}.
#'
#' @details For \code{v_star}, the forest structure doesn't need to be complete.
#' That is,
#' in \code{C}, some trivial intervals \code{c(i,i)}
#' corresponding to regions that are atoms may be missing.
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @examples
#' m <- 20
#' C <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(
#'     c(4, 5), c(8, 9), c(10, 10), c(12, 12),
#'     c(13, 15), c(17, 17), c(18, 19)
#'   ),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' ZL <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#' v_star(1, C, ZL, leaf_list)
#'
#' v_star(1:5, C, ZL, leaf_list)
#'
#' v_star(13:15, C, ZL, leaf_list)
#'
#' v_star(1:m, C, ZL, leaf_list)
#' @export
v_star <- function(S, C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  vec <- numeric(nb_leaves)
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in C and ZL
  # this initialization also takes care of the minima
  # between \zeta_k and card(S inter R_k)
  for (i in 1:nb_leaves) {
    vec[i] <- sum(leaf_list[[i]][length(leaf_list[[i]])] >= S)
  }
  vec <- vec - c(0, vec[1:(nb_leaves - 1)])
  # this way of initializing is much faster than the following naive approach,
  # the downside is that it heavily uses
  # that the elements of leaf_list are sorted
  # naive approach: for (i in 1:nb_leaves) vec[i] <- sum(S %in% leaf_list[[i]])
  for (h in H:1) {
    nb_regions <- length(C[[h]])
    if (nb_regions > 0) {
      for (k in 1:nb_regions) {
        C_HK <- C[[h]][[k]]
        sum_succ <- sum(vec[C_HK[1]:C_HK[2]])
        res <- min(ZL[[h]][k], sum_succ)
        vec[C_HK[1]:C_HK[2]] <- 0
        vec[C_HK[1]] <- res
      }
    }
  }
  return(sum(vec))
}

#' Prune a forest structure to speed up computations
#'
#' @description The pruned forest structure makes the computation of [v_star()],
#' [curve_v_star_forest_naive()] and [curve_v_star_forest_fast()] faster.
#'
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#' @param ZL A list of integer vectors representing the upper
#' bounds \eqn{\zeta_k} of the forest structure.
#' See [v_star()] for more information.
#' @param leaf_list A list of vectors representing the atoms
#' of the forest structure. See [v_star()] for more information.
#' @param prune_leafs A boolean, \code{FALSE} by default.
#' If \code{TRUE}, will also prune atoms/leafs for
#' which \eqn{\zeta_k \geq |R_k|},
#' this makes the computation of [v_star()] and [curve_v_star_forest_naive()]
#' even faster but should not be used with [curve_v_star_forest_fast()]
#' because this needs the structure to be complete (i.e., with all its atoms).
#' This is why the default option is \code{FALSE}.
#' @param delete_gaps A boolean, \code{FALSE} by default.
#' If \code{TRUE}, will also delete the gaps in the structure induced
#' by the pruning, see [delete_gaps()].
#'
#' @return A list with three named elements. \describe{
#' \item{\code{VstarNm}}{\eqn{V^*(\mathbb N_m)} is computed as by-product by
#' the algorithm, so we might as well return it.}
#' \item{\code{C}}{The new \code{C} after pruning.}
#' \item{\code{ZL}}{The new \code{ZL} after pruning.}
#' }
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). Fast confidence bounds for the false discovery
#' proportion over a path of hypotheses.
#' Computo.
#' @export
pruning <- function(C, ZL,
                    leaf_list,
                    prune_leafs = FALSE,
                    delete_gaps = FALSE) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  vec <- numeric(nb_leaves)
  for (i in 1:nb_leaves) {
    vec[i] <- length(leaf_list[[i]])
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
        C_HJ <- C[[h]][[j]]
        candidate <- ZL[[h]][j]
        sum_succ <- sum(vec[C_HJ[1]:C_HJ[2]])
        if (candidate >= sum_succ) {
          res <- sum_succ
          if (prune_leafs || (C_HJ[1] < C_HJ[2])) {
            # Either the region is not a leaf,
            # or it's a leaf and we chose to prune the leafs
            # so we prune.
            ZL[[h]] <- ZL[[h]][-j]
            C[[h]][[j]] <- NULL
          } else {
            # The region is a leaf/atom and we chose
            # to keep the leafs, so we keep it,
            # but we change the \zeta_i to the length of P_i.
            ZL[[h]][j] <- sum_succ
          }
        } else {
          res <- candidate
        }
        vec[C_HJ[1]:C_HJ[2]] <- 0
        vec[C_HJ[1]] <- res
      }
    }
  }
  if (delete_gaps) {
    gaps_deleted <- delete_gaps(C, ZL, leaf_list)
    C <- gaps_deleted$C
    ZL <- gaps_deleted$ZL
  }
  return(list(
    VstarNm = sum(vec),
    C = C,
    ZL = ZL
  ))
}

#' Delete the gaps induced by pruning
#'
#' @description
#' A small optimization that can be done after pruning and
#' that can speed up computations
#' (it removes the gaps introduced by the pruning)
#'
#' @details
#' See [pruning()].
#'
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#' @param ZL A list of integer vectors representing the upper
#' bounds \eqn{\zeta_k} of the forest structure.
#' See [v_star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the
#' forest structure. See [v_star()] for more information.
#'
#' @return A list with two named elements. \describe{
#' \item{\code{C}}{The new \code{C} after deleting the gaps.}
#' \item{\code{ZL}}{The new \code{ZL} after deleting the gaps.}
#' }
#'
#' @usage delete_gaps(C, ZL, leaf_list)
#'
#' @export
delete_gaps <- function(C, ZL, leaf_list) {
  H <- length(C)
  nb_leaves <- length(leaf_list)
  continue <- TRUE
  NEW_C <- list()
  NEW_ZL <- list()
  loop_counter <- 1
  while (continue) {
    NEW_C[[loop_counter]] <- list()
    NEW_ZL[[loop_counter]] <- numeric(0)
    leaf_available <- !logical(nb_leaves)
    regions_delete <- list()
    for (h in 1:H) {
      nb_regions <- length(C[[h]])
      if (nb_regions > 0) {
        for (l in 1:nb_regions) {
          C_HL <- C[[h]][[l]]
          if (all(leaf_available[C_HL[1]:C_HL[2]])) {
            NEW_C[[loop_counter]][[C_HL[1]]] <- C_HL
            NEW_ZL[[loop_counter]][C_HL[1]] <- ZL[[h]][[l]]
            leaf_available[C_HL[1]:C_HL[2]] <- FALSE
            regions_delete <- append(regions_delete, list(c(h, l)))
          }
        }
      }
    }
    for (couple in rev(regions_delete)) {
      h <- couple[1]
      l <- couple[2]
      C[[h]][[l]] <- NULL
      ZL[[h]] <- ZL[[h]][-l]
    }
    if (length(NEW_C[[loop_counter]]) > 0) {
      for (l in rev(seq_along(NEW_C[[loop_counter]]))) {
        if (is.null(NEW_C[[loop_counter]][[l]])) {
          NEW_C[[loop_counter]][[l]] <- NULL
          NEW_ZL[[loop_counter]] <- NEW_ZL[[loop_counter]][-l]
        }
      }
    }
    loop_counter <- loop_counter + 1
    continue <- any(sapply(X = ZL, FUN = function(vec) {
      length(vec) > 0
    }))
  }
  return(list(C = NEW_C, ZL = NEW_ZL))
}

#' Compute a curve of post hoc bounds based on a reference
#' family with forest structure
#'
#' @name curve_v_star_forest
#'
#' @description
#' Computes the post hoc upper bound \eqn{V^*(S_t)} on the number
#' of false positives in a
#' given sequence of selection sets \eqn{S_t} of hypotheses, such that
#' \eqn{S_t\subset S_{t+1}} and \eqn{|S_t| = t},
#' using a reference family \eqn{(R_k, \zeta_k)}
#' that possess the forest structure
#' (see References).
#'
#' @details Two functions are available
#' \describe{
#' \item{\code{curve_v_star_forest_naive}}{
#' Repeatedly calls [v_star()] on each \eqn{S_t},
#' which is not optimized and time-consuming,
#' this should not be used in practice.
#' }
#' \item{\code{curve_v_star_forest_fast}}{
#' A fast and optimized version that leverage the fact that\eqn{S_{t+1}} is
#' the union of \eqn{S_t} and a single hypothesis index.
#' The algorithm needs to work on a complete forest,
#' so this version first completes the forest
#' (unless told that the forest has already been completed,
#' see [forest_completion()]),
#' and the completion fails if the input is a pruned forest (see [pruning()]),
#' so if a pruned forest is given as input,
#' it MUST be said with the \code{is_pruned} argument
#' so that the function skips completion
#' (so the pruned forest given as input must also be complete).
#' }
#' }
#'
#' @param perm An integer vector of elements in \code{1:m}, all different,
#' and of size up to \code{m}
#' (in which case it's a permutation, hence the name).
#' The set \eqn{S_t} is represented by \code{perm[1:t]}.
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds
#' \eqn{\zeta_k} of the forest structure. See [v_star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the
#' forest structure. See [v_star()] for more information.
#' @param pruning A boolean, \code{FALSE} by default.
#' Whether to prune the forest (see [pruning()]) before computing the bounds.
#' Ignored if \code{is_pruned} is \code{TRUE}.
#' @param is_pruned A boolean, \code{FALSE} by default. If \code{TRUE},
#' assumes that the forest structure has already been completed
#' (see [forest_completion()]) and then pruned (see [pruning()])
#' and so skips the completion step and optional pruning step.
#' Must be set to \code{TRUE} if giving a pruned forest, see Details.
#' @param is_complete A boolean, \code{FALSE} by default.
#' If \code{TRUE}, assumes that the forest structure has already been completed
#' (see [forest_completion()]) and so skips the completion step.
#' Ignored if \code{is_pruned} is \code{TRUE}.
#' @param delete_gaps A boolean, \code{FALSE} by default.
#' If \code{TRUE}, will also delete the gaps in the structure induced
#' by the pruning, see [delete_gaps()].
#' Ignored if \code{pruning} is \code{FALSE}.
#'
#' @return A vector of length of same length as \code{perm},
#' where the \code{t}-th
#' element is \eqn{V^*(S_t)}.
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). Fast confidence bounds for the
#' false discovery proportion over a path of hypotheses.
#' Computo.
#' @export
#' @examples
#' m <- 20
#' C <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(
#'     c(4, 5), c(8, 9), c(10, 10), c(12, 12),
#'     c(13, 15), c(17, 17), c(18, 19)
#'   ),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' ZL <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#' curve_v_star_forest_naive(1:m, C, ZL, leaf_list, pruning = FALSE)
#'
#' curve_v_star_forest_naive(1:m, C, ZL, leaf_list, pruning = TRUE)
#'
#' curve_v_star_forest_fast(1:m, C, ZL, leaf_list, pruning = FALSE)
#'
#' curve_v_star_forest_fast(1:m, C, ZL, leaf_list, pruning = TRUE)

#' @rdname curve_v_star_forest
#' @export
curve_v_star_forest_naive <- function(perm, C, ZL,
                                      leaf_list,
                                      pruning = FALSE,
                                      delete_gaps = FALSE) {
  vstars <- numeric(length(perm))

  # the naive version doesn't need a proper completion of the
  # forest structure because v_star
  # implicitly completes, and for the same reason
  # it is allowed to prune leaves

  if (pruning) {
    pruned <- pruning(C, ZL, leaf_list,
      prune_leafs = TRUE,
      delete_gaps = delete_gaps
    )
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

  S <- numeric(0)
  j <- 0
  for (t in perm) {
    j <- j + 1
    S <- c(S, t)
    vstars[j] <- v_star(S, C, ZL, leaf_list)
  }
  return(vstars)
}

#' @rdname curve_v_star_forest
#' @export
curve_v_star_forest_fast <- function(perm, C, ZL,
                                     leaf_list,
                                     pruning = FALSE,
                                     is_pruned = FALSE,
                                     is_complete = FALSE,
                                     delete_gaps = FALSE) {
  vstars <- numeric(length(perm))

  if (!is_pruned) {
    if (!is_complete) {
      # the fast version needs a proper completion of the
      # forest structure, and for the same reason
      # it must not prune the leaves
      completed <- forest_completion(C, ZL, leaf_list)
      C <- completed$C
      ZL <- completed$ZL
    }

    if (pruning) {
      is_pruned <- TRUE # useless atm, but kept for clarity
      pruned <- pruning(C, ZL, leaf_list,
        prune_leafs = FALSE, delete_gaps = delete_gaps
      )
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
  K_MINUS <- list()
  M <- matrix(0, ncol = H, nrow = m)
  for (h in 1:H) {
    zeta_depth_h <- ZL[[h]]
    length_zeta_depth_h <- length(zeta_depth_h)
    etas[[h]] <- rep(0, length_zeta_depth_h)
    K_MINUS[[h]] <- vector("list", length(C[[h]]))
    if (length_zeta_depth_h > 0) {
      for (k in 1:length_zeta_depth_h) {
        if (zeta_depth_h[k] == 0) {
          K_MINUS[[h]][[k]] <- C[[h]][[k]]
        }
        first_leaf <- leaf_list[[C[[h]][[k]][1]]]
        last_leaf <- leaf_list[[C[[h]][[k]][2]]]
        M[first_leaf[1]:last_leaf[length(last_leaf)], h] <- k
      }
    }
  }

  for (t in seq_along(perm)) {
    i_t <- perm[t]
    if (t > 1) {
      previous_vstar <- vstars[t - 1]
    } else {
      previous_vstar <- 0
    }

    # SEARCHING IF i_t IS IN K_MINUS
    # if so, we let go_next = TRUE
    # and we just go next to step t+1
    go_next <- FALSE
    for (h in 1:H) {
      k <- M[i_t, h]
      if ((k > 0) && (!is.null(K_MINUS[[h]][[k]]))) {
        go_next <- TRUE
        break
      }
    }

    # COMPUTING VSTAR AND UPDATING K_MINUS AND ETAS
    if (go_next) {
      # Here, i_t is in K_MINUS
      vstars[t] <- previous_vstar
    } else {
      # Here, i_t isn't in K_MINUS
      for (h in 1:H) {
        k <- M[i_t, h]
        if (k > 0) {
          # if k == 0,
          # there is no k^{(t,h)} because either there is a
          # gap in the structure (maybe because of pruning, or because
          # of the way the structure has been built),
          # either because we are past the max depth of i_t;
          # in this case we just don't do anything
          etas[[h]][[k]] <- etas[[h]][[k]] + 1
          if (etas[[h]][[k]] >= ZL[[h]][[k]]) {
            K_MINUS[[h]][[k]] <- C[[h]][[k]]
            break
          }
        }
      }
      vstars[t] <- previous_vstar + 1
    }
  }
  return(vstars)
}

#' Complete a forest structure
#'
#' @description Completes the forest in the sens of the Reference:
#' adds the missing atoms/leafs
#' in the reference family with a forest structure \eqn{(R_k, \zeta_k)}
#' so that each atom is well represented by a \eqn{R_k}.
#' The associated \eqn{\zeta_k} is
#' taken as the trivial \eqn{|R_k|}.
#'
#' @details The forest must not be pruned (with [pruning()]) beforehand.
#' The code will not behave as expected
#' and will return a wrong result if a pruned forest is given as input.
#' Maybe the function could be rewritten
#' going from the leaves to the roots instead of the contrary,
#' to avoid this issue.
#'
#' @param C A list of list representing the forest structure.
#' See [v_star()] for more information.
#' @param ZL A list of integer vectors representing the upper bounds
#' \eqn{\zeta_k} of the forest structure. See [v_star()] for more information.
#' @param leaf_list A list of vectors representing the atoms of
#' the forest structure. See [v_star()] for more information.
#'
#' @return A list with two named elements. \describe{
#' \item{\code{C}}{The new \code{C} after completion.}
#' \item{\code{ZL}}{The new \code{ZL} after completion.}
#' }
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). A fast algorithm to compute a curve of
#' confidence upper bounds for the False Discovery Proportion using a reference
#' family with a forest structure. arXiv:2502.03849.
#' @export
forest_completion <- function(C, ZL, leaf_list) {
  H <- length(C)

  leaves_to_place <- seq_along(leaf_list)
  len_to_place <- length(leaves_to_place)
  to_delete <- numeric(0)

  for (h in 1:H) {
    j <- 1
    l <- 1

    while (j <= len_to_place) {
      expected_leaf <- leaves_to_place[j]
      end_of_line <- l > length(C[[h]])
      if (!end_of_line) {
        C_HL <- C[[h]][[l]]
      }
      if ((!end_of_line) && expected_leaf == C_HL[[1]]) {
        if (expected_leaf == C_HL[[2]]) {
          to_delete <- c(to_delete, j)
        }
        j <- j + C_HL[2] - C_HL[1] + 1
      } else {
        C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)), l - 1)
        ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]), l - 1)
        to_delete <- c(to_delete, j)
        j <- j + 1
      }
      l <- l + 1
    }

    if (length(to_delete) > 0) {
      leaves_to_place <- leaves_to_place[-to_delete]
      len_to_place <- length(leaves_to_place)
      to_delete <- numeric(0)
    }
  }

  if (len_to_place > 0) {
    h <- H + 1
    C[[h]] <- list()
    ZL[[h]] <- numeric(0)
    for (expected_leaf in leaves_to_place) {
      C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)))
      ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]))
    }
  }

  return(list(C = C, ZL = ZL))
}

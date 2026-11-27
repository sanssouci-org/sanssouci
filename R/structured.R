#' Estimate the number of true null hypotheses
#' among a deterministic set of p-values
#'
#' @description
#' An upper bound on the number of true null hypotheses
#' in a deterministic region
#' is computed with confidence \code{1 - lambda},
#' using the the \eqn{p}-values \code{p_values} associated to
#' those hypotheses.
#' The functions described here can be used as the \code{method} argument
#' of [structured_zetas()].
#'
#' @param p_values A vector of \eqn{p}-values.
#' @param lambda A numeric value in \eqn{(0,1)},
#' the target error level of the test.
#' @param k The positive integer \eqn{k} used by the
#' \eqn{k}-Bonferroni procedure. By default it is equal to 1.
#' @param ... Additional arguments that may be passed
#' to specific \code{zeta} functions.
#' @name structured_zeta
#' @return The number of true nulls is over-estimated as follows:
#' \describe{
#' \item{\code{structured_zeta_dkwm}}{
#' Inversion of the Dvoretzky-Kiefer-Wolfowitz-Massart inequality
#' with parameter \code{lambda}
#' }
#' \item{\code{structured_zeta_hb}}{
#' Number of conserved hypotheses by the Holm-Bonferroni
#' procedure with parameter \code{lambda}
#' }
#' \item{\code{structured_zeta_kbonf}}{
#' Number of conserved hypotheses by the \eqn{k}-Bonferroni
#' procedure with parameter \code{lambda}, plus \eqn{k-1}
#' }
#' \item{\code{structured_zeta_trivial}}{
#' The size of the p-value set which is the trivial
#' upper bound (\code{lambda} is not used in that case)
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
#' @references Durand G. (2025). Fast confidence bounds for the false
#' discovery proportion over a path of hypotheses. Computo.
#' @examples
#' x <- rnorm(100, mean = c(rep(c(0, 2), each = 50)))
#' p_values <- 1 - pnorm(x)
#' lambda <- 0.05
#' structured_zeta_trivial(p_values, lambda)
#'
#' structured_zeta_hb(p_values, lambda)
#'
#' structured_zeta_dkwm(p_values, lambda)
#'
#' structured_zeta_kbonf(p_values, lambda)
#'
#' structured_zeta_kbonf(p_values, lambda, k = 2)
NULL

#' @rdname structured_zeta
#' @export
structured_zeta_hb <- function(p_values, lambda, ...) {
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

#' @rdname structured_zeta
#' @export
structured_zeta_kbonf <- function(p_values, lambda, k = 1, ...) {
  m <- length(p_values)
  threshold <- lambda * k / m
  indexes <- which(p_values <= threshold)
  bound <- m - length(indexes) + (k - 1)
  return(min(bound, m))
}

#' @rdname structured_zeta
#' @export
structured_zeta_trivial <- function(p_values, lambda, ...) {
  return(length(p_values))
}

#' @rdname structured_zeta
#' @export
structured_zeta_dkwm <- function(p_values, lambda, ...) {
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
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#'
#' @return An integer, the number of regions.
#'
#' @export
n_elements_forest_structure <- function(forest_structure) {
  H <- length(forest_structure)
  count <- 0
  for (h in H:1) {
    count <- count + length(forest_structure[[h]])
  }
  return(count)
}

#' Estimate of the proportion of true nulls in each node of a forest structure
#'
#' @description
#' Takes a forest structure as input, given by
#' the couple \code{forest_structure}, \code{leaf_list}
#' and returns the corresponding \eqn{\zeta_k}'s
#' according to the method(s) chosen.
#'
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the forest
#' structure. See [pred_subset_structured()] for more information.
#' @param method A function with arguments \code{(p_values, lambda, ...)}
#' that can
#' compute an upper bound on the false positives in the region associated to
#' the \eqn{p}-values \code{p_values} at confidence level \code{1 - lambda}.
#' It can also be a list of such functions, where the \code{h}-th function
#' is used at depth \code{h} in the tree structure, that is to the \eqn{R_k}'s
#' represented by the elements found in \code{forest_structure[[h]]}.
#' Finally, it can also be a list
#' of lists of such functions, mimicking the structure
#' of \code{forest_structure} itself,
#' that is, \code{method[[h]][[j]]} is applied to the \eqn{R_k}
#' represented by \code{forest_structure[[h]][[j]]}.
#' @param p_values A vector of \eqn{p}-values, must be of size \code{m},
#' with \code{m} the highest element found in the vectors of \code{leaf_list}.
#' @param alpha A target error level in \eqn{(0,1)}.
#' @param stepdown A boolean, \code{FALSE} by default.
#' Whether to use the step-down refinement to try to
#' produce smaller \eqn{\zeta_k}'s, see Details.
#' @param verbose A boolean, \code{FALSE} by default.
#' Whether to print information about the (possibly multiple)
#' round(s) of step-down refinement.
#' @param ... Additional arguments that may be passed to
#' specific \code{structured_zeta} functions.
#'
#' @return \code{zetas_list}: A list of integer vectors representing
#' the upper bounds \eqn{\zeta_k} of the forest structure.
#' See [pred_subset_structured()] for more information.
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
#' forest_structure <- dd$forest_structure
#' method <- structured_zeta_trivial
#' zetas_list <- structured_zetas(
#'   forest_structure,
#'   leaf_list,
#'   method,
#'   p_values,
#'   alpha = 0.05
#' )
#' zetas_list
structured_zetas <- function(forest_structure,
                             leaf_list,
                             method,
                             p_values,
                             alpha,
                             stepdown = FALSE,
                             verbose = FALSE,
                             ...) {
  H <- length(forest_structure)
  K <- n_elements_forest_structure(forest_structure)
  zetas_list <- list()
  NEW_K <- K
  continue <- TRUE
  n_loop <- 0
  while (continue) {
    USAGE_K <- NEW_K
    NEW_K <- K
    for (h in H:1) {
      current_depth <- forest_structure[[h]]
      len <- length(current_depth)
      zeta_inter <- numeric(len)
      if (len > 0) {
        for (k in 1:len) {
          current_node <- current_depth[[k]]
          extracted_p_values <- p_values[unlist(
            leaf_list[current_node[1]:current_node[2]]
          )]
          args_zeta <- list(...)
          if ("p_cdfs" %in% names(args_zeta)) {
            args_zeta$p_cdfs <- args_zeta$p_cdfs[unlist(
              leaf_list[current_node[1]:current_node[2]]
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
      zetas_list[[h]] <- zeta_inter
    }
    if (verbose) {
      n_loop <- n_loop + 1
      print(paste0(
        "loop number=",
        n_loop,
        ", USAGE_K=",
        USAGE_K,
        ", NEW_K=",
        NEW_K
      ))
    }
    continue <- stepdown && (NEW_K < USAGE_K) && (NEW_K > 0)
  }
  return(zetas_list)
}

#' False positives upper bound on a single subset
#' for a structured reference family
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
#' @param forest_structure A list of list representing the forest structure.
#' Each list of \code{forest_structure} represents a level of
#' depth in the forest structure.
#' So, \code{forest_structure[[1]]} lists the regions at depth 1,
#' \code{forest_structure[[2]]} lists the regions at depth 2,
#' and so on. We then use the fact that each region of the reference family
#' is the union of some of its atom over an integer interval.
#' So he elements of each list \code{forest_structure[[i]]}
#' are integer vectors of length 2
#' representing this interval.
#' For example, \code{forest_structure[[1]][[1]] = c(1, 5)}
#' means that the first region
#' at depth 1 is the union of the first 5 atoms,
#' \code{forest_structure[[2]][[3]] = c(4, 5)} means
#' that the third region at depth 2
#' is the union of the atoms 4 and 5, and
#' \code{forest_structure[[3]][[5]] = c(5, 5)}
#' means that the fifth region at depth 3
#' is simply the fifth atom.
#' @param zetas_list A list of integer vectors representing the upper bound
#' \eqn{\zeta_k} associated to a region \eqn{R_k} in the reference family.
#' \code{zetas_list[[h]][j]} is the \eqn{\zeta_k} associated to the
#' \eqn{R_k} described by \code{forest_structure[[h]][[j]]}.
#' @param leaf_list A list of vectors. Each vector is an integer array.
#' The i-th vector contains the indices
#' of the hypotheses in the i-th atom.
#' Atoms form a partition of the set of hypotheses indices :
#' there cannot be overlap, and each index has to be inside one of the atoms.
#' It is really important, for all functions
#' using \code{leaf_list} in this package, that each vector is sorted,
#' and furthermore that
#' \code{leaf_list[[i]]==leaf_list[[i]][1]:leaf_list[[i]][2]}, that is,
#' the leaves follow a 1-dimensional structure.
#' In particular, some functions of the package
#' implicitly use that
#' \code{leaf_list[[i]][1]} is the lowest hypothesis
#' index in leaf \eqn{P_i} and that
#' \code{leaf_list[[i]][length(leaf_list[[i]])]} is the largest.
#'
#' @return An integer value, the post hoc upper bound \eqn{V^*(S)}.
#'
#' @details For \code{pred_subset_structured},
#' the forest structure doesn't need to be complete.
#' That is,
#' in \code{forest_structure},
#' some intervals \code{c(i,i)}
#' corresponding to regions that are atoms may be missing.
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @examples
#' m <- 20
#' forest_structure <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(
#'     c(4, 5), c(8, 9), c(10, 10), c(12, 12),
#'     c(13, 15), c(17, 17), c(18, 19)
#'   ),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' zetas_list <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#' pred_subset_structured(1, forest_structure, zetas_list, leaf_list)
#'
#' pred_subset_structured(1:5, forest_structure, zetas_list, leaf_list)
#'
#' pred_subset_structured(13:15, forest_structure, zetas_list, leaf_list)
#'
#' pred_subset_structured(1:m, forest_structure, zetas_list, leaf_list)
#' @export
pred_subset_structured <- function(S, forest_structure, zetas_list, leaf_list) {
  H <- length(forest_structure)
  n_leaves <- length(leaf_list)
  vec <- numeric(n_leaves)
  # the initialization term for each atom P_i
  # is equivalent to completing the family if it isn't,
  # assuming that leaf_list does indeed contain all leaves
  # and some were just eventually missing in forest_structure and zetas_list
  for (i in 1:n_leaves) {
    vec[i] <- sum(S <= leaf_list[[i]][length(leaf_list[[i]])])
  }
  vec <- vec - c(0, vec[1:(n_leaves - 1)])
  # this way of initializing is much faster than the following naive approach,
  # the downside is that it heavily uses
  # that the elements of leaf_list are sorted
  # naive approach: for (i in 1:n_leaves) vec[i] <- sum(S %in% leaf_list[[i]])
  for (h in H:1) {
    # starting at H takes care of the minima
    # between zeta_k and card(S inter R_k)
    for (k in seq_along(forest_structure[[h]])) {
      current_node <- forest_structure[[h]][[k]]
      sum_succ <- sum(vec[current_node[1]:current_node[2]])
      res <- min(zetas_list[[h]][k], sum_succ)
      vec[current_node[1]:current_node[2]] <- 0
      vec[current_node[1]] <- res
    }
  }
  return(sum(vec))
}

#' Prune a forest structure to speed up computations
#'
#' @description The pruned forest structure makes the computation
#' of [pred_subset_structured()] and
#' [pred_pathway_structured()] faster.
#'
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param zetas_list A list of integer vectors representing the upper
#' bounds \eqn{\zeta_k} of the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param leaf_list A list of vectors representing the atoms/leaves
#' of the forest structure. See [pred_subset_structured()] for more information.
#' @param prune_leafs A boolean, \code{FALSE} by default.
#' If \code{TRUE}, will also prune atoms/leaves for
#' which \eqn{\zeta_k \geq |R_k|},
#' this makes the computation of [pred_subset_structured()]
#' even faster but should not be used with [pred_pathway_structured()]
#' because this needs the structure to be complete (i.e., with all its atoms).
#' This is why the default option is \code{FALSE}.
#' @param delete_gaps A boolean, \code{FALSE} by default.
#' If \code{TRUE}, will also delete the gaps in the structure induced
#' by the pruning, see [delete_gaps()].
#'
#' @return A list with three named elements. \describe{
#' \item{\code{v_star_all}}{\eqn{V^*(\mathbb N_m)} is computed as a by-product
#' by the algorithm, so we might as well return it.}
#' \item{\code{forest_structure}}{The new \code{forest_structure}
#' after pruning.}
#' \item{\code{zetas_list}}{The new \code{zetas_list} after pruning.}
#' }
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). Fast confidence bounds for the false discovery
#' proportion over a path of hypotheses.
#' Computo.
#' @export
pruning <- function(forest_structure,
                    zetas_list,
                    leaf_list,
                    prune_leafs = FALSE,
                    delete_gaps = FALSE) {
  H <- length(forest_structure)
  n_leaves <- length(leaf_list)
  vec <- numeric(n_leaves)
  for (i in 1:n_leaves) {
    vec[i] <- length(leaf_list[[i]])
    # The initialization term for each atom/leaf P_i
    # is equivalent to completing the family if it isn't,
    # assuming that leaf_list does indeed contain all leaves
    # and some were just eventually missing in forest_structure and zetas_list.
    # Using the length of the atoms ensures that we
    # will catch atoms P_i with a \zeta_i larger than its length
  }
  for (h in H:1) {
    for (j in rev(seq_along(forest_structure[[h]]))) {
      current_node <- forest_structure[[h]][[j]]
      current_zeta <- zetas_list[[h]][j]
      sum_succ <- sum(vec[current_node[1]:current_node[2]])
      if (current_zeta >= sum_succ) {
        res <- sum_succ
        if (prune_leafs || (current_node[1] < current_node[2])) {
          # Either the region is not a leaf,
          # or it's a leaf and we chose to prune the leaves
          # so we prune.
          zetas_list[[h]] <- zetas_list[[h]][-j]
          forest_structure[[h]][[j]] <- NULL
        } else {
          # The region is a leaf/atom and we chose
          # to keep the leaves, so we keep it,
          # but we change the \zeta_i to be the length of P_i.
          zetas_list[[h]][j] <- sum_succ
        }
      } else {
        res <- current_zeta
      }
      vec[current_node[1]:current_node[2]] <- 0
      vec[current_node[1]] <- res
    }
  }
  if (delete_gaps) {
    gaps_deleted <- delete_gaps(forest_structure, zetas_list, leaf_list)
    forest_structure <- gaps_deleted$forest_structure
    zetas_list <- gaps_deleted$zetas_list
  }
  return(list(
    v_star_all = sum(vec),
    forest_structure = forest_structure,
    zetas_list = zetas_list
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
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param zetas_list A list of integer vectors representing the upper
#' bounds \eqn{\zeta_k} of the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the
#' forest structure. See [pred_subset_structured()] for more information.
#'
#' @return A list with two named elements. \describe{
#' \item{\code{forest_structure}}{The new \code{forest_structure}
#' after deleting the gaps.}
#' \item{\code{zetas_list}}{The new \code{zetas_list} after deleting the gaps.}
#' }
#'
#' @usage delete_gaps(forest_structure, zetas_list, leaf_list)
#'
#' @export
delete_gaps <- function(forest_structure, zetas_list, leaf_list) {
  H <- length(forest_structure)
  n_leaves <- length(leaf_list)
  continue <- TRUE
  new_forest_structure <- list()
  new_zetas_list <- list()
  loop_counter <- 1
  while (continue) {
    new_forest_structure[[loop_counter]] <- list()
    new_zetas_list[[loop_counter]] <- numeric(0)
    leaf_available <- !logical(n_leaves)
    regions_delete <- list()
    for (h in 1:H) {
      for (l in seq_along(forest_structure[[h]])) {
        current_node <- forest_structure[[h]][[l]]
        if (all(leaf_available[current_node[1]:current_node[2]])) {
          new_forest_structure[[loop_counter]][[current_node[1]]] <-
            current_node
          new_zetas_list[[loop_counter]][current_node[1]] <-
            zetas_list[[h]][[l]]
          leaf_available[current_node[1]:current_node[2]] <- FALSE
          regions_delete <- append(regions_delete, list(c(h, l)))
        }
      }
    }
    for (couple in rev(regions_delete)) {
      h <- couple[1]
      l <- couple[2]
      forest_structure[[h]][[l]] <- NULL
      zetas_list[[h]] <- zetas_list[[h]][-l]
    }
    # TODO: instead of adding elements in new_forest_structure at index
    # current_node[1], add them at index idx_counter and then increment
    # idx_counter, so that the following "if" becomes useless
    # and can be safely removed?
    if (length(new_forest_structure[[loop_counter]]) > 0) {
      for (l in rev(seq_along(new_forest_structure[[loop_counter]]))) {
        if (is.null(new_forest_structure[[loop_counter]][[l]])) {
          new_forest_structure[[loop_counter]][[l]] <- NULL
          new_zetas_list[[loop_counter]] <- new_zetas_list[[loop_counter]][-l]
        }
      }
    }
    loop_counter <- loop_counter + 1
    continue <- any(sapply(X = zetas_list, FUN = function(vec) {
      length(vec) > 0
    }))
  }
  return(list(
    forest_structure = new_forest_structure,
    zetas_list = new_zetas_list
  ))
}

#' Compute a curve of post hoc bounds based on a reference
#' family with forest structure
#'
#' @name pred_pathway_structured
#'
#' @description
#' Computes the post hoc upper bound \eqn{V^*(S_t)} on the number
#' of false positives in a
#' given sequence of selection sets \eqn{S_t} of hypotheses, such that
#' \eqn{S_t\subset S_{t+1}} and \eqn{|S_t| = t}
#' (we call such a sequence a pathway),
#' using a reference family \eqn{(R_k, \zeta_k)}
#' that possess the forest structure
#' (see References).
#'
#' @details A fast and optimized algorithm that leverages
#' the fact that\eqn{S_{t+1}} is
#' the union of \eqn{S_t} and of a single hypothesis index.
#' The algorithm needs to work on a complete forest,
#' so this implementation first completes the forest
#' (unless told that the forest has already been completed,
#' see [forest_completion()]),
#' and the completion fails if the input is a pruned forest (see [pruning()]),
#' so if a pruned forest is given as input,
#' it MUST be said so with the \code{is_pruned} argument
#' so that the function skips completion
#' (hence the pruned forest given as input must already be complete).
#'
#' @param perm An integer vector of elements in \code{1:m}, all different,
#' and of size up to \code{m}
#' (in which case it's a permutation, hence the name).
#' The set \eqn{S_t} is represented by \code{perm[1:t]}.
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param zetas_list A list of integer vectors representing the upper bounds
#' \eqn{\zeta_k} of the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param leaf_list A list of vectors representing the atoms of the
#' forest structure. See [pred_subset_structured()] for more information.
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
#' @return A vector of length equal to the length of \code{perm},
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
#' forest_structure <- list(
#'   list(c(2, 5), c(8, 15), c(16, 19)),
#'   list(c(3, 5), c(8, 10), c(12, 15), c(16, 16), c(17, 19)),
#'   list(
#'     c(4, 5), c(8, 9), c(10, 10), c(12, 12),
#'     c(13, 15), c(17, 17), c(18, 19)
#'   ),
#'   list(c(8, 8), c(9, 9), c(13, 13), c(14, 15), c(18, 18), c(19, 19))
#' )
#' zetas_list <- list(
#'   c(4, 8, 4),
#'   c(3, 3, 4, 1, 3),
#'   c(2, 2, 1, 1, 2, 1, 2),
#'   c(1, 1, 1, 2, 1, 1)
#' )
#' leaf_list <- as.list(1:m)
#'
#' pred_pathway_structured(1:m, forest_structure,
#' zetas_list, leaf_list, pruning = FALSE)
#'
#' pred_pathway_structured(1:m, forest_structure,
#' zetas_list, leaf_list, pruning = TRUE)

#' @rdname pred_pathway_structured
#' @export
pred_pathway_structured <- function(
  perm,
  forest_structure,
  zetas_list,
  leaf_list,
  pruning = FALSE,
  is_pruned = FALSE,
  is_complete = FALSE,
  delete_gaps = FALSE
) {
  vstars <- numeric(length(perm))

  if (!is_pruned) {
    if (!is_complete) {
      # the fast version needs a proper completion of the
      # forest structure, and for the same reason
      # it must not prune the leaves
      completed <- forest_completion(forest_structure, zetas_list, leaf_list)
      forest_structure <- completed$forest_structure
      zetas_list <- completed$zetas_list
    }

    if (pruning) {
      is_pruned <- TRUE # useless atm, but kept for clarity
      pruned <- pruning(
        forest_structure,
        zetas_list,
        leaf_list,
        prune_leafs = FALSE,
        delete_gaps = delete_gaps
      )
      forest_structure <- pruned$forest_structure
      zetas_list <- pruned$zetas_list
      m <- length(unlist(leaf_list))

      if (length(perm) == m) {
        # means that the last bound to compute
        # is V^*({1, ..., m}),
        # but the pruning already computed
        # V^*({1, ..., m}) as a by-product so we
        # might as well use it:
        vstars[m] <- pruned$v_star_all
        perm <- perm[-m]
      }
    }
  }

  H <- length(forest_structure)
  m <- length(unlist(leaf_list))

  # preparation of the initial value of the etas (a copy of zetas but
  # with only zeroes), of the initial value of K^- (with the R_k's
  # such that zeta_k = 0) and of a m x H matrix M such that
  # M[i, h] gives the index k of forest_structure[[h]] such that hypothesis i
  # is in the R_k given by forest_structure[[h]][[k]]
  etas <- list()
  K_MINUS <- list()
  M <- matrix(0, nrow = m, ncol = H)
  for (h in 1:H) {
    zeta_depth_h <- zetas_list[[h]]
    etas[[h]] <- rep(0, length(zeta_depth_h))
    K_MINUS[[h]] <- vector("list", length(forest_structure[[h]]))
    for (k in seq_along(zeta_depth_h)) {
      if (zeta_depth_h[k] == 0) {
        K_MINUS[[h]][[k]] <- forest_structure[[h]][[k]]
      }
      first_leaf <- leaf_list[[forest_structure[[h]][[k]][1]]]
      last_leaf <- leaf_list[[forest_structure[[h]][[k]][2]]]
      M[first_leaf[1]:last_leaf[length(last_leaf)], h] <- k
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
          if (etas[[h]][[k]] >= zetas_list[[h]][[k]]) {
            K_MINUS[[h]][[k]] <- forest_structure[[h]][[k]]
            break
          }
        }
      }
      vstars[t] <- previous_vstar + 1
    }
  }
  return(vstars)
}

# TODO: impose by design and by outputs that all forest_structures
# are without gaps (no more optional delete_gaps argument), in which
# case forest_completion should work after pruning?
# even if we prune the leafs???

#' Complete a forest structure
#'
#' @description Completes the forest in the sense of the Reference:
#' adds the missing atoms/leaves
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
#' @param forest_structure A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param zetas_list A list of integer vectors representing the upper bounds
#' \eqn{\zeta_k} of the forest structure.
#' See [pred_subset_structured()] for more information.
#' @param leaf_list A list of vectors representing the atoms of
#' the forest structure. See [pred_subset_structured()] for more information.
#'
#' @return A list with two named elements. \describe{
#' \item{\code{forest_structure}}{The new \code{forest_structure}
#' after completion.}
#' \item{\code{zetas_list}}{The new \code{zetas_list} after completion.}
#' }
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025). A fast algorithm to compute a curve of
#' confidence upper bounds for the False Discovery Proportion using a reference
#' family with a forest structure. arXiv:2502.03849.
#' @export
forest_completion <- function(forest_structure, zetas_list, leaf_list) {
  H <- length(forest_structure)

  leaves_to_place <- seq_along(leaf_list)
  len_to_place <- length(leaves_to_place)
  to_delete <- numeric(0)

  for (h in 1:H) {
    j <- 1
    l <- 1

    while (j <= len_to_place) {
      expected_leaf <- leaves_to_place[j]
      end_of_line <- l > length(forest_structure[[h]])
      if (!end_of_line) {
        current_node <- forest_structure[[h]][[l]]
      }
      if ((!end_of_line) && expected_leaf == current_node[[1]]) {
        if (expected_leaf == current_node[[2]]) {
          to_delete <- c(to_delete, j)
        }
        j <- j + current_node[2] - current_node[1] + 1
      } else {
        forest_structure[[h]] <- append(
          forest_structure[[h]],
          list(c(expected_leaf, expected_leaf)),
          after = l - 1
        )
        zetas_list[[h]] <- append(
          zetas_list[[h]],
          length(leaf_list[[expected_leaf]]),
          after = l - 1
        )
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
    forest_structure[[h]] <- list()
    zetas_list[[h]] <- numeric(0)
    for (expected_leaf in leaves_to_place) {
      forest_structure[[h]] <- append(
        forest_structure[[h]],
        list(c(expected_leaf, expected_leaf))
      )
      zetas_list[[h]] <- append(
        zetas_list[[h]],
        length(leaf_list[[expected_leaf]])
      )
    }
  }

  return(list(forest_structure = forest_structure, zetas_list = zetas_list))
}

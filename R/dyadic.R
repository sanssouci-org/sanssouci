#' Create a complete dyadic tree structure
#'
#' @description
#' Produce a set of regions with forest structure,
#' with a single tree, which is dyadic.
#'
#'
#' @name dyadic
#' @details \describe{
#' \item{\code{dyadic_from_leaf_list}}{
#' Dyadic tree structure from a given list of atoms/leafs
#' }
#' \item{\code{dyadic_from_window_size}}{
#' Dyadic tree structure from window size:
#' the number of elements in each atom/leaf is set to \code{s}
#' }
#' \item{\code{dyadic_from_height}}{
#' Dyadic tree structure from height:
#' the total height of the tree is set to \code{H}
#' }
#' }
#'
#' @param flavor An integer equal to 1 or 2 (which is the default).
#' If \code{flavor == 1}, start from the leaves
#' and group nodes of a same height 2 by 2 as
#' long as possible (note that this can introduce gaps).
#' If \code{flavor==2}, start from the root and divide nodes in
#' 2 nodes of equal
#' size as long as possible.
#'
#' @return A list with two named elements:\describe{
#' \item{\code{leaf_list}}{
#' A list of vectors representing the atoms/leaves of the forest structure.
#' See [pred_subset_structured()] for more information.
#' }
#' \item{\code{forest_structure}}{
#' A list of list representing the forest structure.
#' See [pred_subset_structured()] for more information.
#' }
#' }
#'
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Durand G. (2025).
#' Fast confidence bounds for the false discovery
#' proportion over a path of hypotheses. Computo.
#' @examples
#' m <- 6
#' dd <- dyadic_from_window_size(m, s = 2, flavor = 2)
#' str(dd)
#'
#' dd <- dyadic_from_height(m, H = 3, flavor = 2)
#' str(dd)
#'
#' dd <- dyadic_from_height(m, flavor = 2)
#' str(dd)
#'
#' leaf_list <- dd$leaf_list
#' dd <- dyadic_from_leaf_list(leaf_list, flavor = 2)
#' str(dd)
NULL

#' @param leaf_list A list of vectors representing the atoms of the
#' forest structure. See [pred_subset_structured()] for more information.
#' @export
#' @rdname dyadic
dyadic_from_leaf_list <- function(leaf_list, flavor = 2) {
  nb_leaves <- length(leaf_list)
  if (flavor == 1) {
    inter <- seq_len(nb_leaves)
    current_depth <- mapply(c, inter, inter, SIMPLIFY = FALSE)
    forest_structure <- list(current_depth)
    repeat {
      len <- length(current_depth)
      if (len == 1) {
        break
      }
      new_depth <- list()
      j <- 1
      while (j < len) {
        new_depth <- c(
          new_depth,
          list(c(current_depth[[j]][1], current_depth[[j + 1]][2]))
        )
        j <- j + 2
      }
      current_depth <- new_depth
      forest_structure <- c(list(current_depth), forest_structure)
    }
  } else if (flavor == 2) {
    current_depth <- list(c(1, nb_leaves))
    forest_structure <- list(current_depth)
    continue <- (nb_leaves > 1)
    while (continue) {
      continue <- FALSE
      old_depth <- current_depth
      current_depth <- list()
      len <- length(old_depth)
      for (i in seq_len(len)) {
        node <- old_depth[[i]]
        leaves_in_node <- node[2] - node[1] + 1
        if (leaves_in_node > 1) {
          cut_in_2 <- ceiling(leaves_in_node / 2)
          current_depth <- c(
            current_depth,
            list(c(node[1], node[1] + cut_in_2 - 1)),
            list(c(node[1] + cut_in_2, node[2]))
          )
          if (cut_in_2 > 1) {
            continue <- TRUE
          }
        }
      }
      forest_structure <- c(forest_structure, list(current_depth))
    }
  }
  return(list(leaf_list = leaf_list, forest_structure = forest_structure))
}

#' @param m An integer value, the number of hypotheses to have in the structure
#' @param s An integer value, the number of elements in each atom/leaf.
#' Default value is 1.
#' @export
#' @rdname dyadic
dyadic_from_window_size <- function(m, s = 1, flavor = 2) {
  leaves <- floor(m / s)
  leaf_list <- list()
  for (l in 1:leaves) {
    leaf <- seq_len(s) + (l - 1) * s
    leaf_list <- c(leaf_list, list(leaf))
  }
  if (s * leaves < m) {
    leaf <- seq(1 + l * s, m)
    leaf_list <- c(leaf_list, list(leaf))
  }
  forest_structure <- dyadic_from_leaf_list(leaf_list, flavor)$forest_structure
  return(list(leaf_list = leaf_list, forest_structure = forest_structure))
}

#' @param H An integer value, the desired maximal height of the tree.
#' If NULL (by default),
#' use \code{floor(2 + log2(m - 1))} (or 1 if \code{m} = 1), which is the
#' maximum achievable height given \code{m}.
#' @export
#' @rdname dyadic
dyadic_from_height <- function(
  m,
  H = ifelse(m == 1, 1, floor(2 + log2(m - 1))),
  flavor = 2
) {
  if (m <= 2^(H - 2)) {
    OLD_H <- H
    H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
    warning(
      "H=",
      OLD_H,
      " is too large for m=",
      m,
      ", \nH=",
      OLD_H,
      " reduced to H=",
      H
    )
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
      leaf_list <- c(
        leaf_list,
        list(plus1 * (base + 1) + seq_len(base) + (i - 1) * base)
      )
    }
  }
  forest_structure <- dyadic_from_leaf_list(leaf_list, flavor)$forest_structure
  return(list(leaf_list = leaf_list, forest_structure = forest_structure))
}

# the forest must not be pruned beforehand
# the following code fails if the input is a pruned forest
# this version is really cool because 
# it doesn't add atoms 1 and 2 if R_(1,2) is here
# and both atoms are missing, which is indeed useless
# if we apply our V.star computing functions
# buuuut the pruning function may then delete 
# R_(1,2) and produce an incomplete forest,
# that is because the pruning function is only 
# built to protect real atoms
# so we have to keep this one hidden forever,
# or until we find a use for it,
# and use a version that really add all atoms
# on the bright side, that makes the code simpler
# TODO BEFORE MERGE: move to a "useless function atm" file?
forest.completion.hidden <- function(C, ZL, leaf_list) {
  H <- length(C)
  
  leaves.to.place <- 1:length(leaf_list)
  len.to.place <- length(leaves.to.place)
  to.delete <- numeric(0)
  
  for (h in 1:H) {
    
    j <- 1
    l <- 1
    l_ancestor <- 1
    
    while (j <= len.to.place) {
      expected_leaf <- leaves.to.place[j]
      end_of_line <- l > length(C[[h]])
      if (! end_of_line) {Chl <- C[[h]][[l]]}
      while ((! end_of_line) && h > 1 && Chl[1] > C[[h-1]][[l_ancestor]][2]) {
        l_ancestor <- l_ancestor + 1
      }
      
      if ((! end_of_line) && expected_leaf == Chl[[1]]) {
        if (expected_leaf == Chl[[2]]) {
          to.delete <- c(to.delete, j)
        }
        j <- j + Chl[2] - Chl[1] +1
        l <- l + 1
      }
      else {
        if (h == 1 || C[[h-1]][[l_ancestor]][1] < Chl[1]) {
          C[[h]] <- append(C[[h]], list(c(expected_leaf, expected_leaf)), l - 1)
          ZL[[h]] <- append(ZL[[h]], length(leaf_list[[expected_leaf]]), l - 1)
          l <- l + 1
        }
        to.delete <- c(to.delete, j)
        j <- j + 1
      }
      
    }
    
    leaves.to.place <- leaves.to.place[-to.delete]
    len.to.place <- length(leaves.to.place)
    to.delete <- numeric(0)
  }
  
  
  return(list(C = C, ZL = ZL))
}

# useless at the moment because finally
# curve.V.star.forest.fast does not need K.1
# TODO BEFORE MERGE: move to a "useless function atm" file?
# adapt this code and apply it recursively to delete gaps in a completed
# and pruned forest? sounds like it should work
compute.K.1 <- function(C, ZL, leaf_list) {
  H <- length(C)
  K.1 <- list()
  nb_leaves <- length(leaf_list)
  leaf.available <- ! logical(nb_leaves)
  for (h in 1:H) {
    K.1[[h]] <- list()
    nb_regions <- length(C[[h]])
    if (nb_regions > 0) {
      for (l in 1:nb_regions) {
        Chl <- C[[h]][[l]]
        if (all(leaf.available[Chl[1]:Chl[2]])) {
          K.1[[h]][[l]] <- Chl
          leaf.available[Chl[1]:Chl[2]] <- FALSE
        }
      }
    }
  }
  return(K.1)
}

# TODO BEFORE MERGE: move to a "useless function atm" file?
nodeLabel <- function(x) {
  paste(unlist(x), collapse = ":")
}

# TODO later: understand and maybe revamp?
# boolean: is the set x included in y?
setinclude <- function(x, y) {
  all(is.element(x, y))
}

# TODO later: understand and maybe revamp?
# internal function of tree.from.list computes a 'proto' tree structure from a list of pairs region/zeta
proto.tree.from.list <- function(listR) {
  lengths <- sapply(listR, function(lll) {
    length(lll[[1]])
  })
  o <- order(lengths, decreasing = TRUE)
  listR <- listR[o]
  protoCZ <- list(list(listR[[1]]))
  # protoZ<-list()
  lenlist <- length(o)
  parkourlist <- 2
  while (parkourlist <= lenlist) {
    curr_elem <- listR[[parkourlist]]
    b1 <- FALSE
    hb <- 0
    lb <- 0
    H <- length(protoCZ)
    for (h in H:1) {
      if (b1) {
        break
      }
      for (l in 1:length(protoCZ[[h]])) {
        if (setinclude(curr_elem$R, protoCZ[[h]][[l]]$R)) {
          b1 <- TRUE
          hb <- h
          lb <- l
          break
        }
      }
    }
    if (b1) {
      if (hb + 1 <= H) {
        b2 <- TRUE
        for (p in 1:length(protoCZ[[hb + 1]])) {
          for (q in 1:length(protoCZ[[hb]])) {
            if (setinclude(protoCZ[[hb + 1]][[p]]$R, protoCZ[[hb]][[q]]$R)) 
              break
          }
          if (q >= lb) {
            protoCZ[[hb + 1]] <- append(protoCZ[[hb + 1]], list(curr_elem), p - 1)
            b2 <- FALSE
            break
          }
        }
        if (b2) 
          protoCZ[[hb + 1]] <- c(protoCZ[[hb + 1]], list(curr_elem))
      } else {
        protoCZ <- c(protoCZ, list(list(curr_elem)))
      }
    } else {
      protoCZ[[1]] <- c(protoCZ[[1]], list(curr_elem))
    }
    parkourlist <- parkourlist + 1
  }
  return(protoCZ)
}

# TODO later: understand and maybe revamp?
# internal function of tree.from.list recursive function designed to compute a list of leaves
recurleaves <- function(numvect, trunc) {
  leaf_list <- list()
  is.end <- (length(trunc) == 1)
  start <- TRUE
  Union <- numeric(0)
  for (j in 1:length(trunc[[1]])) {
    curr_vect <- trunc[[1]][[j]]$R
    # recherche des mecs de trunc[[1]] qui sont dans numvect
    if (setinclude(curr_vect, numvect)) {
      if (start) {
        start <- FALSE
      } # useless if, could just be start <- FALSE
      Union <- union(Union, curr_vect)
      if (is.end) {
        leaf_list <- c(leaf_list, list(curr_vect))
      } else {
        leaf_list <- c(leaf_list, recurleaves(curr_vect, trunc[-1]))
      }
    } else {
      if (!start) 
        break
    }
  }
  if (start) {
    leaf_list <- list(numvect)
  } else {
    compl <- setdiff(numvect, Union)
    if (length(compl) > 0) {
      leaf_list <- c(leaf_list, list(compl))
    }
  }
  return(leaf_list)
}

# TODO later: understand and maybe revamp?
# computes an incomplete tree structure, a tree of zetas, and a list of leaves from a list of pairs region/zeta.
tree.from.list <- function(m, listR) {
  protoCZ <- proto.tree.from.list(listR)
  leaf_list <- recurleaves(1:m, protoCZ)
  leaves <- length(leaf_list)
  H <- length(protoCZ)
  C <- vector("list", length = H)
  ZL <- vector("list", length = H)
  for (h in 1:H) {
    parkourleaves <- 1
    leaf_start <- 1
    len <- length(protoCZ[[h]])
    C[[h]] <- vector("list", length = len)
    for (j in 1:len) {
      start <- TRUE
      ZL[[h]] <- c(ZL[[h]], protoCZ[[h]][[j]]$z)
      curr_elem <- protoCZ[[h]][[j]]$R
      while ((parkourleaves <= leaves) && (start || setinclude(leaf_list[[parkourleaves]], curr_elem))) {
        if (start && setinclude(leaf_list[[parkourleaves]], curr_elem)) {
          leaf_start <- parkourleaves
          start <- FALSE
        }
        parkourleaves <- parkourleaves + 1
      }
      C[[h]][[j]] <- c(leaf_start, parkourleaves - 1)
    }
  }
  return(list(C = C, ZL = ZL, leaf_list = leaf_list))
}

# makes a complete dydadic tree structure C from a list of leaves method 1 starts from the leaves and groups
# nodes of a same height 2 by 2 as long as possible method 2 starts from the root and divide nodes in 2 equal
# nodes as long as possible
dyadic.from.leaf_list <- function(leaf_list, method) {
    leafs <- length(leaf_list)
    if (method == 1) {
        inter <- seq_len(leafs)
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
        Ch <- list(c(1, leafs))
        C <- list(Ch)
        continue <- TRUE
        while (continue) {
            continue <- FALSE
            oldCh <- Ch
            Ch <- list()
            len <- length(oldCh)
            for (i in seq_len(len)) {
                leafs_in_node <- oldCh[[i]][2] - oldCh[[i]][1] + 1
                if (leafs_in_node > 1) {
                  cut2 <- ceiling(leafs_in_node/2)
                  Ch <- c(Ch, list(c(oldCh[[i]][1], oldCh[[i]][1] + cut2 - 1)), list(c(oldCh[[i]][1] + cut2, oldCh[[i]][2])))
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
# makes a complete dyadic tree structure and a list of leaves from the desired size s of the leaves
dyadic.from.window.size <- function(m, s, method) {
    leafs <- floor(m/s)
    leaf_list <- list()
    for (l in 1:leafs) {
        leaf <- seq_len(s) + (l - 1) * s
        leaf_list <- c(leaf_list, list(leaf))
    }
    if (s * leafs < m) {
        leaf <- seq(1 + l * s, m)
        leaf_list <- c(leaf_list, list(leaf))
        # leafs<-leafs+1
    }
    C <- dyadic.from.leaf_list(leaf_list, method)
    return(list(leaf_list = leaf_list, C = C))
}
# makes a complete dyadic tree structure and a list of leaves from the desired maximal height H
dyadic.from.height <- function(m, H, method) {
    if (m <= 2^(H - 2)) {
        oldH <- H
        H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
        warning("H=", oldH, " is too large for m=", m, ", \nH=", oldH, " reduced to H=", H)
    }
    if (m < 2^(H - 1)) {
        leaf_list <- as.list(1:m)
    } else {
        leaf_list <- list()
        base <- m%/%2^(H - 1)
        plus1 <- m%%2^(H - 1)
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
# makes a complete dyadic tree structure and a list of leaves from the maximum achievable height H
dyadic.max.height <- function(m, method) {
    H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
    return(dyadic.from.height(m, H, method))
}
# computes zeta from an HB test at level lambda
zeta.HB <- function(pval, lambda) {
    s <- length(pval)
    k <- 0
    sp <- sort(pval)
    CONT <- TRUE
    while ((k < s) && CONT) {
        if (sp[k + 1] > lambda/(s - k)) {
            CONT <- FALSE
        } else {
            k <- k + 1
        }
    }
    # res<-list(hatk = k, rej_ind = sort(or[seq_len(k)]) )
    return(s - k)
}
# computes zeta equals to the size of the set
zeta.trivial <- function(pval, lambda) {
    return(length(pval))
}
# computes zeta from the DKWM/Storey method, at level lambda
zeta.DKWM <- function(pval, lambda) {
    s <- length(pval)
    sorted.pval <- c(0, sort(pval))
    dkwm <- min((sqrt(log(1/lambda)/2)/(2 * (1 - sorted.pval)) + sqrt(log(1/lambda)/(8 * (1 - sorted.pval)^2) + 
        (s - seq(0, s))/(1 - sorted.pval)))^2)
    return(min(s, floor(dkwm)))
}
# number of unique regions in a given tree
nb.elements <- function(C) {
    H <- length(C)
    count <- length(C[[H]])
    if (H > 1) {
        for (h in (H - 1):1) {
            for (j in 1:length(C[[h]])) {
                if (C[[h]][[j]][1] < C[[h]][[j]][2]) 
                  count <- count + 1
            }
        }
    }
    return(count)
}
# computes the tree of the zetas from given tree structure, list of leaves, method to compute the zetas, and
# level alpha, using the union bound method, that is the method to compute the zetas is applied at level lambda
# = alpha / number of regions
zetas.tree <- function(C, leaf_list, method, pvalues, alpha) {
    H <- length(C)
    K <- nb.elements(C)
    leafs <- length(leaf_list)
    zeta_leafs <- numeric(leafs)
    for (i in 1:length(C[[H]])) {
        if (C[[H]][[i]][1] == C[[H]][[i]][2]) {
            zeta_leafs[C[[H]][[i]][1]] <- method(pvalues[leaf_list[[C[[H]][[i]][1]]]], alpha/K)
        }
    }
    ZL <- list()
    for (h in H:1) {
        len <- length(C[[h]])
        zeta_inter <- numeric(len)
        for (j in 1:len) {
            if (C[[h]][[j]][1] < C[[h]][[j]][2]) {
                zeta_inter[j] <- method(pvalues[unlist(leaf_list[C[[h]][[j]][1]:C[[h]][[j]][2]])], alpha/K)
            } else {
                zeta_inter[j] <- zeta_leafs[C[[h]][[j]][1]]
            }
        }
        ZL[[h]] <- zeta_inter
    }
    return(ZL)
}
# computes the tree of the zetas from given tree structure, list of leaves, method to compute the zetas, and
# level alpha, using the refined method, that is the method which tries to estimate the number of regions with
# only signal, if any
zetas.tree.refined <- function(C, leaf_list, method, pvalues, alpha) {
    H <- length(C)
    K <- nb.elements(C)
    leafs <- length(leaf_list)
    zeta_leafs <- numeric(leafs)
    continue <- TRUE
    new_K <- K
    while (continue) {
        usage_K <- new_K
        new_K <- K
        for (i in 1:length(C[[H]])) {
            if (C[[H]][[i]][1] == C[[H]][[i]][2]) {
                zeta_leafs[C[[H]][[i]][1]] <- method(pvalues[leaf_list[[C[[H]][[i]][1]]]], alpha/usage_K)
                if (zeta_leafs[C[[H]][[i]][1]] == 0) 
                  new_K <- new_K - 1
            }
        }
        ZL <- list()
        for (h in H:1) {
            len <- length(C[[h]])
            zeta_inter <- numeric(len)
            for (j in 1:len) {
                if (C[[h]][[j]][1] < C[[h]][[j]][2]) {
                  zeta_inter[j] <- method(pvalues[unlist(leaf_list[C[[h]][[j]][1]:C[[h]][[j]][2]])], alpha/usage_K)
                  if (zeta_inter[j] == 0) 
                    new_K <- new_K - 1
                } else {
                  zeta_inter[j] <- zeta_leafs[C[[h]][[j]][1]]
                }
            }
            ZL[[h]] <- zeta_inter
        }
        if (new_K == usage_K) 
            continue <- FALSE
    }
    return(ZL)
}
# boolean: is the set x included in y?
setinclude <- function(x, y) {
    all(is.element(x, y))
}
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
            if (b1) 
                break
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
# internal function of tree.from.list recursive function designed to compute a list of leaves
recurLeafs <- function(numvect, trunc) {
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
            }
            Union <- union(Union, curr_vect)
            if (is.end) {
                leaf_list <- c(leaf_list, list(curr_vect))
            } else {
                leaf_list <- c(leaf_list, recurLeafs(curr_vect, trunc[-1]))
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
# computes an incomplete tree structure, a tree of zetas, and a list of leaves from a list of pairs region/zeta.
tree.from.list <- function(m, listR) {
    protoCZ <- proto.tree.from.list(listR)
    leaf_list <- recurLeafs(1:m, protoCZ)
    leafs <- length(leaf_list)
    H <- length(protoCZ)
    C <- vector("list", length = H)
    ZL <- vector("list", length = H)
    for (h in 1:H) {
        parkourleafs <- 1
        leaf_start <- 1
        len <- length(protoCZ[[h]])
        C[[h]] <- vector("list", length = len)
        for (j in 1:len) {
            start <- TRUE
            ZL[[h]] <- c(ZL[[h]], protoCZ[[h]][[j]]$z)
            curr_elem <- protoCZ[[h]][[j]]$R
            while ((parkourleafs <= leafs) && (start || setinclude(leaf_list[[parkourleafs]], curr_elem))) {
                if (start && setinclude(leaf_list[[parkourleafs]], curr_elem)) {
                  leaf_start <- parkourleafs
                  start <- FALSE
                }
                parkourleafs <- parkourleafs + 1
            }
            C[[h]][[j]] <- c(leaf_start, parkourleafs - 1)
        }
    }
    return(list(C = C, ZL = ZL, leaf_list = leaf_list))
}
# computes V star (S) from a complete tree structure, a complete zeta tree, and a list of leaves
V.star.all.leafs <- function(S, C, ZL, leaf_list) {
    H <- length(C)
    leafs <- length(leaf_list)
    Vec <- numeric(leafs)
    id <- seq_len(leafs)
    for (i in 1:leafs) {
        Vec[i] <- min(ZL[[H]][i], length(intersect(S, leaf_list[[i]])))
    }
    if (H > 1) {
        for (h in (H - 1):1) {
            len <- length(C[[h]])
            vec_inter <- numeric(leafs)
            new_id <- numeric(len)
            for (j in 1:len) {
                intra <- min(min(ZL[[h]][j], length(intersect(S, unlist(leaf_list[C[[h]][[j]][1]:C[[h]][[j]][2]])))), 
                  sum(Vec[id[which((C[[h]][[j]][1] <= id) & (C[[h]][[j]][2] >= id))]]))
                vec_inter[C[[h]][[j]][1]] <- intra
                new_id[j] <- C[[h]][[j]][1]
            }
            Vec <- vec_inter
            id <- new_id
        }
    }
    return(sum(Vec[id]))
}
# completes an incomplete tree structure
tree.expand <- function(C, ZL, leaf_list) {
    H <- length(C)
    leafs <- length(leaf_list)
    for (i in 1:H) {
        expected_leaf <- 1
        len <- length(C[[i]])
        j <- 1
        while (j <= len) {
            if (expected_leaf == C[[i]][[j]][1]) {
                # TOUT VA BIEN suivant
                expected_leaf <- C[[i]][[j]][2] + 1
                j <- j + 1
            } else {
                # PROBLÈME résolution
                C[[i]] <- append(C[[i]], list(c(expected_leaf, expected_leaf)), j - 1)
                ZL[[i]] <- append(ZL[[i]], length(leaf_list[[expected_leaf]]), j - 1)
                # suivant
                expected_leaf <- expected_leaf + 1
                j <- j + 1
                len <- len + 1
            }
        }
        while (expected_leaf < (leafs + 1)) {
            # PROBLÈME résolution
            C[[i]] <- c(C[[i]], list(c(expected_leaf, expected_leaf)))
            ZL[[i]] <- c(ZL[[i]], length(leaf_list[[expected_leaf]]))
            # suivant
            expected_leaf <- expected_leaf + 1
        }
    }
    if (len == leafs) {
        return(list(C = C, ZL = ZL))
    } else {
        CH <- list()
        ZLH <- numeric(0)
        for (j in 1:len) {
            id_start <- C[[H]][[j]][1]
            id_end <- C[[H]][[j]][2]
            if (id_start == id_end) {
                CH <- c(CH, list(c(id_start, id_start)))
                ZLH <- c(ZLH, ZL[[H]][j])
            } else {
                for (k in id_start:id_end) {
                  CH <- c(CH, list(c(k, k)))
                  ZLH <- c(ZLH, length(leaf_list[[k]]))
                }
            }
        }
        C <- c(C, list(CH))
        ZL <- c(ZL, list(ZLH))
    }
    return(list(C = C, ZL = ZL))
}
# computes V star (S) from an incomplete tree structure, an incomplete zeta tree, and a list of leaves
V.star <- function(S, C, ZL, leaf_list) {
    all_leafs <- tree.expand(C, ZL, leaf_list)
    return(V.star.all.leafs(S, all_leafs$C, all_leafs$ZL, leaf_list))
}

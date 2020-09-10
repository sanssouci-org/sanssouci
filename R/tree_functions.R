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
#' dd <- dyadic.from.max.height(m, method=2)
#' str(dd)
#'
#' leaf_list <- dd$leaf_list

NULL


#' Create a complete dyadic tree structure
#' 
#' @param leaf_list A list of leaves
#' @param method A numeric value. If \code{method == 1}, start from the leaves
#'   and group nodes of a same height 2 by 2 as long as possible. If
#'   \code{method==2}, start from the root and divide nodes in 2 nodes of equal
#'   size as long as possible
#' @return A list of lists containing the dyadic structure
#' 
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

#' @inheritParams dyadic.from.leaf_list
#' @param m An integer value, the number of elements in the structure
#' @param s An integer value, the number of elements in each leaf
#' @return A list with two elements:\describe{
#' \item{leaf_list}{A list of leaves}
#' \item{C}{A list of lists containing the dyadic structure}
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

#' @inheritParams dyadic.from.window.size 
#' @param H An integer value, the desired maximal height of the tree
#' @export
#' @rdname dyadic

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

#' @inheritParams dyadic.from.height
#' @export
#' @rdname dyadic

dyadic.from.max.height <- function(m, method) {
    H <- ifelse(m == 1, 1, floor(2 + log2(m - 1)))
    return(dyadic.from.height(m, H, method))
}

#' Estimate the number of true null hypotheses among a set of p-values
#'
#' @param pval A vector of \eqn{p}-values
#' @param lambda A numeric value in \eqn{[0,1]}, the target level of the test 
#' @name zeta
#' @examples
#' x <- rnorm(100, mean = c(rep(c(0, 2), each = 50)))
#' pval <- 1-pnorm(x)
#' zeta.trivial(pval)
#' zeta.HB(pval, 0.05)
#' zeta.DKWM(pval, 0.05)
NULL
#' @return The numer of true nulls is estimated as follows:
#' \describe{
#' \item{\code{zeta.DKWM}}{Dvoretzky-Kiefer-Wolfowitz-Massart inequality (related to the Storey estimator of the proportion of true nulls) with parameter \code{lambda}}
#' \item{\code{zeta.HB}}{Holm-Bonferroni test with parameter \code{lambda}}
#' \item{\code{zeta.trivial}}{the size of the p-value set (\eqn{lambda} is not used)}
#' }
#' @references Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. The Annals of Mathematical Statistics, pages 642-669.
#' @references Holm, S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6 (1979), pp. 65-70.
#' @references Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. The Annals of Probability, pages 1269-1283.
#' @references Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(3):479-498.
#> NULL

#' @export
#' @rdname zeta
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

#' @export
#' @rdname zeta
zeta.trivial <- function(pval, lambda) {
    return(length(pval))
}

#' @export
#' @rdname zeta
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
            Ch <- C[[h]]
            for (j in 1:length(Ch)) {
                Chj <- Ch[[j]]
                if (Chj[1] < Chj[2]) 
                    count <- count + 1
            }
        }
    }
    return(count)
}

#' Estimate of the proportion of true nulls in each node of a tree
#' 
#' @param C Tree structure
#' @param leaf_list A list of leaves
#' @param method A function to compute the estimators
#' @param pvalues A vector of \eqn{p}-values
#' @param alpha A target level
#' @details The proportion of true nulls in each node is estimated by an union bound on the regions. That is, the provided method is applied at level \code{alpha/nR} where \code{nR} is the number of regions.
#' @export
#' @examples
#' 
#' m <- 1000 
#' dd <- dyadic.from.window.size(m, s = 10, method = 2)
#' leaf_list <- dd$leaf_list
#' mu <- gen.mu.leaves(m, K1 = 3, d = 1, grouped = FALSE, "const", barmu = 4, leaf_list)
#' pvalues<-gen.p.values(m, mu, rho = 0)
#' C <- dd$C 
#' ZL<-zetas.tree(C, leaf_list, zeta.DKWM, pvalues, alpha = 0.05)
zetas.tree <- function(C, leaf_list, method, pvalues, alpha) {
    H <- length(C)
    K <- nb.elements(C)
    leaves <- length(leaf_list)
    zeta_leaves <- numeric(leaves)
    CH <- C[[H]]
    for (i in 1:length(CH)) {
        CHi <- CH[[i]]
        if (CHi[1] == CHi[2]) {
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


#' @rdname zetas.tree
#' 
#' @details In \code{zetas.tree.refined}, one tries to estimate the number of regions containing only signal.
#' 
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
                Chj <- Chj
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
            }
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


V.star.all.leaves <- function(S, C, ZL, leaf_list) {
    H <- length(C)
    leaves <- length(leaf_list)
    Vec <- numeric(length = leaves)
    id <- seq_len(length.out = leaves)
    for (i in 1:leaves) {
#        len <- length(intersect(S, leaf_list[[i]]))
        len <- sum(S %in% leaf_list[[i]])
        Vec[i] <- min(ZL[[H]][i], len)
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
                len <- sum(S %in% lvs)
                rm(lvs)
                ss <- sum(Vec[id[which((Chj[1] <= id) & (Chj[2] >= id))]])
                intra <- min(ZL[[h]][j], len, ss)
                vec_inter[Chj[1]] <- intra
                new_id[j] <- Chj[1]
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
    leaves <- length(leaf_list)
    for (i in 1:H) {
        expected_leaf <- 1
        len <- length(C[[i]])
        j <- 1
        while (j <= len) {
            Cij <- C[[i]][[j]]
            if (expected_leaf == Cij[1]) {
                # TOUT VA BIEN suivant
                expected_leaf <- Cij[2] + 1
                j <- j + 1
            } else {
                # PROBLEME resolution
                C[[i]] <- append(C[[i]], list(c(expected_leaf, expected_leaf)), j - 1)
                ZL[[i]] <- append(ZL[[i]], length(leaf_list[[expected_leaf]]), j - 1)
                # suivant
                expected_leaf <- expected_leaf + 1
                j <- j + 1
                len <- len + 1
            }
        }
        while (expected_leaf < (leaves + 1)) {
            # PROBLEME resolution
            C[[i]] <- c(C[[i]], list(c(expected_leaf, expected_leaf)))
            ZL[[i]] <- c(ZL[[i]], length(leaf_list[[expected_leaf]]))
            # suivant
            expected_leaf <- expected_leaf + 1
        }
    }
    if (len == leaves) {
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

#' Post hoc bound on the number of false positives
#' 
#' @param S A subset of hypotheses
#' @param C Tree
#' @param ZL Zeta tree
#' @param leaf_list List of leaves
#' @return An integer value, upper bound on the number false positives in S
#' @export
V.star <- function(S, C, ZL, leaf_list) {
    all_leaves <- tree.expand(C, ZL, leaf_list)
    return(V.star.all.leaves(S, 
                            all_leaves$C, 
                            all_leaves$ZL, 
                            leaf_list))
}

nodeLabel <- function(x) {
    paste(unlist(x), collapse = ":")
}


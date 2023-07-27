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


#' Create a complete dyadic tree structure
#' 
#' @param leaf_list A list of leaves
#' @param method A numeric value. If \code{method == 1}, start from the leaves
#'   and group nodes of a same height 2 by 2 as long as possible. If
#'   \code{method==2}, start from the root and divide nodes in 2 nodes of equal
#'   size as long as possible
#' @return A list of lists containing the dyadic structure
#' @rdname dyadic
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
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @references Dvoretzky, A., Kiefer, J., and Wolfowitz, J. (1956). Asymptotic minimax character of the sample distribution function and of the classical multinomial estimator. The Annals of Mathematical Statistics, pages 642-669.
#' @references Holm, S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6 (1979), pp. 65-70.
#' @references Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz inequality. The Annals of Probability, pages 1269-1283.
#' @references Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(3):479-498.
#> NULL

#' @export
#' @rdname zeta
zeta.HB <- function(pval, lambda) {
	m <- length(pval)
	sorted.pval <- sort(pval)
	
	thresholds <- lambda / (m - m:1 + 1)
	v <- sorted.pval - thresholds
	indexes <- which(v > 0)
	if (! length(indexes)) {
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
																																			(s - seq(0, s))/(1 - sorted.pval)))^2,
							na.rm=TRUE)
	return(min(s, floor(dkwm)))
}

# number of unique regions in a given tree
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

# number of unique regions in a given tree,
# assuming we didn't extend the leaves
# TODO BEFORE MERGE: rename nb.elements
nb.elements.no.extension <- function(C) {
	H <- length(C)
	count <- 0
	for (h in H:1) {
		count <- count + length(C[[h]])
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
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
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

# TODO BEFORE MERGE: rename zetas.tree, change call of nb.elements, complete documentation
# with refinement
zetas.tree.no.extension <- function(C, leaf_list, method, pvalues, alpha, refine=FALSE, verbose=FALSE) {
	H <- length(C)
	K <- nb.elements.no.extension(C)
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
				if (typeof(method) == "list") {
					if (typeof(method[[h]]) == "list") {
						zeta_method <- method[[h]][[j]]
					} else {
						zeta_method <- method[[h]]
					}
				} else {
					zeta_method <- method
				}
				zeta_inter[j] <- zeta_method(pvals, alpha / usage_K)
				if (refine && (zeta_inter[j] == 0) )
					new_K <- new_K - 1
			}
			ZL[[h]] <- zeta_inter
		}
		if (verbose) {
			nb_loop <- nb_loop + 1
			print(paste0("loop number=", nb_loop,", usage_K=",usage_K,", new_K=",new_K))
		}
		continue <- refine && (new_K < usage_K)
	}
	return(ZL)
}


#' @rdname zetas.tree
#' 
#' @details In \code{zetas.tree.refined}, one tries to estimate the number of regions containing only signal.
#'
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

# NOTE: assumes the forest is extended
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

# the following assumes that zeta_k is <= to
# the cardinal of R_k but this should always happen
# super.prune can also prune atoms for which zeta_k is =
# the cardinal of R_k, this can speed up the
# computation of V.star.no.extension
# but must not be used in conjunction with 
# the fast version of curve.V.star.forest,
# hence it is at FALSE by default
# TODO later: fast curve.Vmax is slower with the pruning
# which is unexpected, maybe it is because the pruning leaves gaps
# => maybe revamp to delete gaps?
# not sure if super.prune and delete.gaps can be both TRUE at the moment
pruning <- function(C, ZL, leaf_list, super.prune = FALSE, delete.gaps = FALSE) {
	H <- length(C)
	nb_leaves <- length(leaf_list)
	Vec <- numeric(nb_leaves) 
	for (i in 1:nb_leaves) {
		Vec[i] <- length(leaf_list[[i]])
	}
	# the initialization term for each atom P_i
	# is equivalent to completing the family if it isn't,
	# assuming that leaf_list does indeed contain all leaves
	# and some were just eventually missing in C and ZL
	for (h in H:1) {
		nb_regions <- length(C[[h]])
		if (nb_regions > 0) {
			for (j in nb_regions:1) {
				Chj <- C[[h]][[j]]
				if (Chj[1]==Chj[2]) { 
					res <- ZL[[h]][j]
					if (super.prune && res >= length(leaf_list[[Chj[1]]])) {
						ZL[[h]] <- ZL[[h]][-j]
						C[[h]][[j]] <- NULL
					}
				} else {
					sum_succ <- sum(Vec[Chj[1]:Chj[2]])
					if (ZL[[h]][j] >= sum_succ) {
						res <- sum_succ
						ZL[[h]] <- ZL[[h]][-j]
						C[[h]][[j]] <- NULL
					} else {
						res <- ZL[[h]][j]
					}
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

# TODO BEFORE MERGE: document
delete.gaps <- function(C, ZL, leaf_list) {
	H <- length(C)
	nb_leaves <- length(leaf_list)
	continue <- TRUE
	newC <- list()
	newZL <- list()
	loop.counter <- 1
	while(continue){
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
			ZL[[h]] <- ZL[[h]][-j]
		}
		if(length(newC[[loop.counter]]) > 0) {
			for(l in length(newC[[loop.counter]]):1) {
				if (is.null(newC[[loop.counter]][[l]])){
					newC[[loop.counter]][[l]] <- NULL
					newZL[[loop.counter]] <- 	newZL[[loop.counter]][-l]
				}
			}
		}
		loop.counter <- loop.counter + 1
		continue <- any(sapply(X = ZL, FUN = function(vec){length(vec) > 0}))
	}
	return(list(C = newC, ZL = newZL))
}

# TODO BEFORE MERGE: change call of V.star, document
curve.V.star.forest.naive <- function(perm, C, ZL, leaf_list, pruning = FALSE, delete.gaps = FALSE){
	vstars <- numeric(length(perm))
	
	# the naive version doesn't need a proper completion of the
	# forest structure because V.star.no.extension
	# implicitly completes, and for the same reason
	# it can use super pruning
	
	if (pruning){
		pruned <- pruning(C, ZL, leaf_list, super.prune = TRUE, delete.gaps = delete.gaps)
		C <- pruned$C
		ZL <- pruned$ZL
		m <- length(unlist(leaf_list))
		if(length(perm) == m){
			# means that length(perm) = m,
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

# the forest must not be pruned beforehand
# the following code fails if the input is a pruned forest
# TODO BEFORE MERGE: document
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
compute.K.1 <- function(C, leaf_list) {
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


# the forest must not be pruned beforehand
# the completion fails if the input is a pruned forest
# TODO BEFORE MERGE: document, change the comments just above that are not accurate anymore
# the forest can't be pruned but incomplete, that is the right condition, so
# if is.pruned we assume that is.complete is also TRUE
curve.V.star.forest.fast <- function(perm, C, ZL, leaf_list, is.pruned = FALSE, is.complete = FALSE, pruning = FALSE, delete.gaps = FALSE){
	
	vstars <- numeric(length(perm))
	
	if (! is.pruned) {
		if (! is.complete) {
			# the fast version needs a proper completion of the
			# forest structure, and for the same reason
			# it must not use super pruning
			completed <- forest.completion(C, ZL, leaf_list)
			C <- completed$C
			ZL <- completed$ZL
		}
		
		if (pruning) {
			is.pruned <- TRUE
			pruned <- pruning(C, ZL, leaf_list, super.prune = FALSE, delete.gaps = delete.gaps)
			C <- pruned$C
			ZL <- pruned$ZL
			m <- length(unlist(leaf_list))
			
			if (length(perm) == m) {
				# means that length(perm) = m,
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

# completes an incomplete tree structure
# NOTE: this automatically results in an extended tree too
# NOTE: this expects an input that is already extended too, if not we can
# add leaves with trivial zeta even if they are already here with proper zeta
# but at smaller depth (doesn't change V.star.all.leaves output luckily)
# NOTE: the expected input is strange given that all the leaves are expected
# to be given by leaf_list, they can only be missing in C, and C has to be
# properly annotated wrt leaf_list and the leaves labels
# TODO BEFORE MERGE: delete this
tree.expand <- function(C, ZL, leaf_list) {
	H <- length(C)
	leaves <- length(leaf_list)
	for (i in 1:H) {
		expected_leaf <- 1 # the label of the leaf we expect to find in the forest structure
		len <- length(C[[i]])
		j <- 1
		while (j <= len) {
			Cij <- C[[i]][[j]]
			if (expected_leaf == Cij[1]) {
				# GOOD, on to the next one
				expected_leaf <- Cij[2] + 1
				j <- j + 1
			} else {
				# PROBLEM => resolution by adding the leaf
				C[[i]] <- append(C[[i]], list(c(expected_leaf, expected_leaf)), j - 1)
				ZL[[i]] <- append(ZL[[i]], length(leaf_list[[expected_leaf]]), j - 1)
				# on to the next one
				expected_leaf <- expected_leaf + 1
				j <- j + 1
				len <- len + 1
			}
			# j <- j + 1 should be here, duplicate code
		}
		while (expected_leaf < (leaves + 1)) {
			# PROBLEM => resolution
			C[[i]] <- c(C[[i]], list(c(expected_leaf, expected_leaf)))
			ZL[[i]] <- c(ZL[[i]], length(leaf_list[[expected_leaf]]))
			# on to the next one
			expected_leaf <- expected_leaf + 1
		}
	}
	if (len == leaves) {
		return(list(C = C, ZL = ZL))
	} else {
		# prepares an additional layer made of all the leaves
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
	return(list(C = C, ZL = ZL)) # duplicate code, change the if condition to (len < leaves)
}

#' Post hoc bound on the number of false positives
#' 
#' @param S A subset of hypotheses
#' @param C Tree
#' @param ZL Zeta tree
#' @param leaf_list List of leaves
#' @return An integer value, upper bound on the number false positives in S
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc false positive control for structured hypotheses. Scandinavian Journal of Statistics, 47(4), 1114-1148.
#' @export
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

# TODO BEFORE MERGE: move to a "useless function atm" file?
nodeLabel <- function(x) {
	paste(unlist(x), collapse = ":")
}


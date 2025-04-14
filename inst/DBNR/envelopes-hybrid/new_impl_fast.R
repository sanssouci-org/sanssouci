curve.V.star.forest.fast.14hypcol <- function(perm, C, ZL, leaf_list, pruning = FALSE, is.pruned = FALSE, is.complete = FALSE, delete.gaps = FALSE){
  
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
  # such that zeta_k = 0) and of a H x m matrix M such that
  # M[h, i] gives the index k of C[[h]] such that hypothesis i
  # is in the R_k given by C[[h]][[k]]
  etas <- ZL
  K.minus <- list()
  M <- matrix(0, ncol = m, nrow = H)
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
        M[h, first_leaf[1]:last_leaf[length(last_leaf)]] <- k
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
      # if (go.next) {
      #   break
      # }
      # for (couple in K.minus[[h]]) {
      #   if (! is.null(couple)) {
      #     lower_leaf <- leaf_list[[couple[1]]]
      #     lower_hyp <- lower_leaf[1]
      #     upper_leaf <- leaf_list[[couple[2]]]
      #     upper_hyp <- upper_leaf[length(upper_leaf)]
      #     if ((i.t >= lower_hyp) && (i.t <= upper_hyp)) {
      #       go.next <- TRUE
      #       # print(paste0(i.t, " is in K minus"))
      #       break
      #     }
      #   }
      # }
    	k <- M[h, i.t]
    	if ((k > 0) && (! is.null(K.minus[[h]][[k]]))){
    		go.next <- TRUE
    		break
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
            # there is no k^{(t,h)} because there is a 
            # gap in the structure (because of pruning)
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

curve.V.star.forest.fast.14hyprow <- function(perm, C, ZL, leaf_list, pruning = FALSE, is.pruned = FALSE, is.complete = FALSE, delete.gaps = FALSE){
  
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
    
    ################################
    # SEARCHING IF i_t IS IN K MINUS
    # if so, go.next == TRUE
    # and we just go next to step t+1
    go.next <- FALSE
    for (h in 1:H) {
      # if (go.next) {
      #   break
      # }
      # for (couple in K.minus[[h]]) {
      #   if (! is.null(couple)) {
      #     lower_leaf <- leaf_list[[couple[1]]]
      #     lower_hyp <- lower_leaf[1]
      #     upper_leaf <- leaf_list[[couple[2]]]
      #     upper_hyp <- upper_leaf[length(upper_leaf)]
      #     if ((i.t >= lower_hyp) && (i.t <= upper_hyp)) {
      #       go.next <- TRUE
      #       # print(paste0(i.t, " is in K minus"))
      #       break
      #     }
      #   }
      # }
    	k <- M[i.t, h]
    	if ((k > 0) && (! is.null(K.minus[[h]][[k]]))){
    		go.next <- TRUE
    		break
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
            # there is no k^{(t,h)} because there is a 
            # gap in the structure (because of pruning)
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
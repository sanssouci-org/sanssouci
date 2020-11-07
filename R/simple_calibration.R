t_inv_linear <- function(y, k, m) {
    y * m / k;
}

#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 6
#' pivStat <- get_pivotal_stat_slow(X, categ, B)
get_pivotal_stat_slow <- function(X, categ, B, 
                               testFUN = t.test,
                               t_inv = t_inv_linear,
                               K = nrow(X)) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values for B permutations of the class assignments
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) ## parametric p-values
    for (bb in 1:B) {
        categ_perm <- sample(categ, length(categ))
        s0 <- which(categ_perm == 0)
        s1 <- which(categ_perm == 1)
        for (ii in 1:m) {
            tt <- testFUN(X[ii, s0], X[ii, s1])
            pval0[ii, bb] <- tt$p.value
        }
    }
    
    ## Step 2: (partially) sort each column
    pval0_sorted <- apply(pval0, 2, sort, partial = 1:K)
    pval0_sorted <- pval0_sorted[1:K, , drop = FALSE]
    
    ## Step 3: apply template function
    tkInv <- matrix(nrow = K, ncol = B)
    for (kk in 1:K) {
        tkInv[kk, ] <- t_inv(pval0_sorted[kk, ], kk, m)
    }
    
    ## Step 4: report min for each column
    matrixStats::colMins(tkInv)
}


#' pivStat <- get_pivotal_stat(X, categ, B)
get_pivotal_stat <- function(X, categ, B, 
                                  rowTestFUN = rowWelchTests,
                                  t_inv = t_inv_linear,
                                  K = nrow(X)) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values for B permutations of the class assignments
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) ## parametric p-values
    for (bb in 1:B) {
        categ_perm <- sample(categ, length(categ))
        rwt <- rowTestFUN(X, categ = categ_perm)
        pval0[, bb] <- rwt$p.value
    }

    ## Step 2: (partially) sort each column
    pval0_sorted <- -partialColSortDesc(-pval0, K)

    ## Step 3: apply template function
    tkInv <- matrix(nrow = m, ncol = B)
    for (kk in 1:K) {
        tkInv[kk, ] <- t_inv(pval0_sorted[kk, ], kk, m)
    }
    
    ## Step 4: report min for each column
    matrixStats::colMins(tkInv)
}

#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pivStat <- get_one_pivotal_stat(X, categ)
#' B <- 10
#' pivStat <- replicate(B, get_one_pivotal_stat_slow(X, sample(categ)))
get_one_pivotal_stat_slow <- function(X, categ, 
                                 testFUN = t.test,
                                 t_inv = t_inv_linear,
                                 K = nrow(X)) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values 
    s0 <- which(categ == 0)
    s1 <- which(categ == 1)
    pval0 <- numeric(m)
    for (ii in 1:m) {
        tt <- testFUN(X[ii, s0], X[ii, s1])
        pval0[ii] <- tt$p.value
    }
    
    ## Step 2: partial sorting
    pval0_sorted <- sort(pval0, partial = 1:K)[1:K]
    
    ## Step 3: apply template function
    tkInv <- numeric(K)
    for (kk in 1:K) {
        tkInv[kk] <- t_inv(pval0_sorted[kk], kk, m)
    }
    
    ## Step 4: report min
    min(tkInv)
}

#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pivStat <- get_one_pivotal_stat(X, categ)
#' B <- 10
#' pivStat <- replicate(B, get_one_pivotal_stat(X, sample(categ)))
get_one_pivotal_stat <- function(X, categ, 
                                 rowTestFUN = rowWelchTests,
                                 t_inv = t_inv_linear,
                                 K = nrow(X)) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values 
    rwt <- rowTestFUN(X, categ = categ)
    pval0 <- rwt$p.value
    
    ## Step 2: partial sorting
    pval0_sorted <- sort(pval0, partial = 1:K)[1:K]
    
    ## Step 3: apply template function
    tkInv <- numeric(K)
    for (kk in 1:K) {
        tkInv[kk] <- t_inv(pval0_sorted[kk], kk, m)
    }
    
    ## Step 4: report min
    min(tkInv)
}
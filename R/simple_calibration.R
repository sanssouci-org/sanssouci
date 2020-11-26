t_inv_linear <- function(y, k, m) {
    y * m / k;
}

t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}

#' @describeIn get_pivotal_stat Get pivotal statistic (slow version)
#' @param testFUN A function with the same I/O as \code{t.test}
#' 
#' @importFrom stats t.test
#' @seealso \code{\link[stats]{t.test}}
#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 6
#' pivStat <- sansSouci:::get_pivotal_stat_slow(X, categ, B)
#' 
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


#' Get pivotal statistic
#' 
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param categ An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#' @param B A numeric value, the number of permutations to be performed
#' @param rowTestFUN a testing function with the same I/O as \code{rowWelchTests}
#' @param t_inv A function with the same I/O as \code{t_inv_linear}
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}. Defaults to \eqn{m}.
#' 
#' @examples
#' 
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 60
#' pivStat <- sansSouci:::get_pivotal_stat_fast(X, categ, B)
get_pivotal_stat_fast <- function(X, categ, B, 
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

#' @describeIn get_pivotal_stat Get one pivotal statistic (slow version)
#' @param testFUN A function with the same I/O as \code{t.test}
#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pivStat <- sansSouci:::get_one_pivotal_stat(X, categ)
#' B <- 10
#' pivStat <- replicate(B, sansSouci:::get_one_pivotal_stat_slow(X, sample(categ)))
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

#' @describeIn get_pivotal_stat Get one pivotal statistic
#' @param testFUN A function with the same I/O as \code{t.test}
#' @examples
#'
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' pivStat <- sansSouci:::get_one_pivotal_stat(X, categ)
#' B <- 10
#' pivStat <- replicate(B, sansSouci:::get_one_pivotal_stat(X, sample(categ)))
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

#' Get permutation p-values
#' 
#' Get a matrix of p-values under the null hypothesis obtained by repeated 
#' permutation of class labels
#' 
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param categ An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#'   specifying the column indices of the first and second samples. If not
#'   provided, a one-sample test is performed.
#' @param B A numeric value, the number of permutations to be performed
#' @param rowTestFUN a testing function with the same I/O as \code{rowWelchTests}
#' @return A \eqn{m \times B} matrix whose entry i,j corresponds to \eqn{p_{(i)}(g_j.X)} with notation of the AoS 2020 paper cited below (section 4.5) 
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @examples
#' 
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' p0 <- sansSouci:::get_perm_p(X, categ, B)
get_perm_p <- function(X, categ, B, 
                       rowTestFUN = rowWelchTests) {
    m <- nrow(X)
    
    ## Step 1: calculate $p$-values for B permutations of the class assignments
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) ## parametric p-values
    for (bb in 1:B) {
        categ_perm <- sample(categ, length(categ))
        rwt <- rowTestFUN(X, categ = categ_perm)
        pval0[, bb] <- rwt$p.value
    }
    
    ## Step 2: sort each column
    return(colSort(pval0))
}


#' Get pivotal statistic
#' 
#' 
#' @param p0 A \eqn{m \times B} matrix of null p-values obtained from \eqn{B} 
#' permutations for \eqn{m} hypotheses
#' @param t_inv A function with the same I/O as \code{t_inv_linear}
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}. Defaults to \eqn{m}.
#' @return A vector of length \eqn{B} pivotal statitics, whose j-th entry corresponds to \eqn{\psi(g_j.X)} with notation of the AoS 2020 paper cited below (section 4.5) 
#'
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' 
#' @examples
#' 
#' set.seed(0xBEEF)
#' m <- 5
#' n <- 45
#' X <- matrix(rnorm(m*n), ncol = n, nrow = m)
#' categ <- rbinom(n, 1, 0.4)
#' B <- 10
#' p0 <- sansSouci:::get_perm_p(X, categ, B)
#' pivStat <- sansSouci:::get_pivotal_stat(p0)
get_pivotal_stat <- function(p0,
                             t_inv = t_inv_linear,
                             K = nrow(p0)) {
    m <- nrow(p0)
    B <- ncol(p0)
    ## Step 3: apply template function
    tkInv <- matrix(nrow = K, ncol = B)
    for (kk in 1:K) {
        tkInv[kk, ] <- t_inv(p0[kk, ], kk, m)
    }
    
    ## Step 4: report min for each column
    matrixStats::colMins(tkInv)
}

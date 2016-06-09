##' Pivotal statistic for the balanced threshold family of BNR
##'
##' @param kmaxH0 A \eqn{c} x \eqn{B} matrix of \code{B} Monte-Carlo
##' samples of \code{c} test statistics under the null hypothesis,
##' with each sample (column) sorted decreasingly.
##' @param m An integer value, the total number of hypotheses
##' tested. \code{m} should not be less than \code{c}.
##' @return the pivotal statistic \eqn{s_k} for the balanced threshold family introduced by BNR
##'
##' @export
##'
kFWERPivotalStatistic1 <- function(kmaxH0, m) {
    B <- ncol(kmaxH0)
    rks <- rowRanks(kmaxH0, ties.method="max") ## 'max' is the default in 'matrixStats'
    colMins(rks)/B
}

kFWERPivotalStatistic2 <- function(kmaxH0, kmaxH0C) {
    B <- ncol(kmaxH0)
    m <- nrow(kmaxH0)
    c <- nrow(kmaxH0C)
    stopifnot(c <= m)
    ## NB: the calculation only requires the first c rows of kmaxH0
    R <- matrix(NA_real_, c, B)
    for (bb in 1:B) {
        for (kk in 1:c) {
            ekb <-kmaxH0[kk, bb]
            ekp <- kmaxH0C[kk, ]
            R[kk, bb] <- mean(ekb >= ekp)
        }
    }
    colMins(R)
}


kFWERPivotalStatistic3 <- function(kmaxH0, kmaxH0C) {
    ## avoid storing all R because we only need the min!
    B <- ncol(kmaxH0)
    m <- nrow(kmaxH0)
    c <- nrow(kmaxH0C)
    stopifnot(c <= m)
    ## NB: the calculation only requires the first c rows of kmaxH0
    res <- rep(Inf, B)
    for (kk in 1:c) {
        ekp <- kmaxH0C[kk, ]
        for (bb in 1:B) {
            rb <- res[bb]
            ekb <-kmaxH0[kk, bb]
            rkb <- mean(ekb >= ekp)
            if (rkb < rb) {
                res[bb] <- rkb
            }
        }
    }
    res
}

kFWERPivotalStatistic4 <- function(kmaxH0, kmaxH0C) {
    B <- ncol(kmaxH0)
    m <- nrow(kmaxH0)
    c <- nrow(kmaxH0C)
    stopifnot(c <= m)
    ## NB: the calculation only requires the first c rows of kmaxH0
    res <- rep(Inf, B)
    for (bb in 1:B) {
        for (kk in 1:c) {
            ekb <-kmaxH0[kk, bb]
            ekp <- kmaxH0C[kk, ]
            rkb <- mean(ekb >= ekp)
            if (rkb < res[bb]) {
                res[bb] <- rkb
            }
        }
    }
    res
}

kFWERPivotalStatistic <- function(mat, kMax, C=1:nrow(mat)) {
    c <- length(C)
    c <- min(kMax, c)  # K \vee |C| in th BNR paper

    kmaxH0 <- partialColSortDesc(mat, c);  ## no need to go further than c (??)
    kmaxH0C <- partialColSortDesc(mat[C, ], c);
    minPseudoRanks(kmaxH0, kmaxH0C)
}

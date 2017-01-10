#' @keywords internal
minPseudoRanksR <- function(kmaxH0, kmaxH0C) {
    B <- ncol(kmaxH0)
    m <- nrow(kmaxH0)
    c <- nrow(kmaxH0C)
    stopifnot(c <= m)
    ## NB: the calculation only requires the first c rows of kmaxH0
    R <- matrix(NA_real_, c, B)
    for (bb in 1:B) {
        for (kk in 1:c) {
            ek <-kmaxH0[kk, ]
            ekb <- kmaxH0C[kk, bb]
            R[kk, bb] <- mean(ek >= ekb)
        }
    }
    colMins(R)
}

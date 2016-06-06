kFWERPivotalStatistic <- function(kmaxH0) {
    B <- ncol(kmaxH0)
    rks <- rowRanks(kmaxH0, ties.method="max") ## 'max' is the default in 'matrixStats'
    colMins(rks)/B
}

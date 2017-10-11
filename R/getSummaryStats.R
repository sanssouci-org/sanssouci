#' Convert a matrix of observations into summary statistics
#' 
#' Convert a matrix of observations labelled into categories into summary 
#' statistics for each category
#' 
#' The following statistics are calculated: sums, sums of squares, means, 
#' standard deviations, sample sizes
#' 
#' @param mat A numeric matrix whose rows correspond to variables and columns to
#'   observations
#'   
#' @param categ A vector of \code{ncol(mat)} categories for the observations
#'   
#' @return A list of \code{ncol(mat)} elements containing the above-described
#'   summary statistics for each category
#' 
#' @author Pierre Neuvial
#'   
#' @examples
#' 
#' if (require("multtest")) { 
#'   data(golub, package="multtest")
#'   stats <- getSummaryStats(golub, categ=golub.cl)
#' } else {
#'   mat <- matrix(rnorm(3051*38), ncol=38)
#'   cls <- rep(c(0, 1), times=c(27, 11))
#'   stats <- getSummaryStats(mat, categ=cls)
#' } 

getSummaryStats <- function(mat, categ=colnames(mat)) {
    stopifnot(ncol(mat) == length(categ))
    cats <- unique(categ)
    res <- list()
    for (cc in seq(along=cats)) {
        ww <- which(categ==cats[cc])
        matc <- mat[, ww]
        
        sumc <- rowSums(matc)
        sum2c <- rowSums(matc^2)
        nc <- length(ww)
        mc <- sumc/nc
        sc <- sqrt((sum2c - sumc^2/nc)/(nc-1))
        
        res[[cc]] <- list(sum = sumc,
                          sum2 = sum2c,
                          n = nc,
                          mean = mc,
                          sd = sc)
    }
    names(res) <- cats
    return(res)
}

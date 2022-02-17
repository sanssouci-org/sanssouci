#' Wilcoxon rank sum tests for each row of a matrix
#'
#' @param mat A \code{m x n} numeric matrix whose rows correspond to variables and columns to
#'   observations
#'
#' @param categ Either a numeric vector of \code{n} categories in \eqn{0, 1} for
#'   the observations, or a \code{n x B} matrix stacking \code{B} such vectors
#'   (typically permutations of an original vector of size \code{n})
#'
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{wilcox.test}}, alternative = "greater" is the alternative that
#'   class 1 is shifted to the right of class 0.
#'
#' @param correct A logical indicating whether to apply continuity correction in
#'   the normal approximation for the p-value.
#'
#' @return A list containing the following components:
#'   \describe{ \item{statistic}{the value of the statistics} \item{p.value}{the
#'   p-values for the tests}} 
#'
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @seealso wilcox.test
#' 
#' @return A list containing the following components:
#' \describe{ 
#'   \item{statistic}{the value of the statistics}
#'   \item{p.value}{the p-values for the tests} 
#'   \item{estimate}{the median difference between groups (only calculated if \code{B=1} for computational efficiency)}}
#'   Each of these elements is a matrix of size \code{m x B}, coerced to a vector of length \code{m} if \code{B=1}
#'   
#' @importFrom matrixStats rowRanks rowTabulates
#' 
#' @details This function performs \code{m x B} Wilcoxon T tests on
#'   \code{n} observations. It is vectorized along the rows of \code{mat}. This
#'   makes the code much faster than using loops of 'apply' functions,
#'   especially for high-dimensional problems (small n and large m) because the
#'   overhead of the call to the 'wilcox.test' function is avoided. Note that it
#'   is not vectorized along the columns of \code{categ} (if any), as a basic
#'   'for' loop is used.
#'   
#' @details The p-values are computed using the normal approximation as
#'   described in the \code{\link{wilcox.test}} function. The exact p-values
#'   (which can be useful for small samples with no ties) are not implemented
#'   (yet).
#'
#' @details For simplicity, 'estimate' returns the difference between the group medians, which does not match the component 'estimate' output by \code{wilcox.test}
#'   

#'
#' @export
#' @examples
#' 
#' p <- 200
#' n <- 50
#' mat <- matrix(rnorm(p*n), ncol = n)
#' cls <- rep(c(0, 1), each = n/2)
#' 
#' stats <- rowWilcoxonTests(mat, categ = cls, alternative = "two.sided")
#' str(stats)
#' 
#' # permutation of class labels
#' cls_perm <- replicate(11, sample(cls))
#' stats <- rowWilcoxonTests(mat, categ = cls_perm, alternative = "two.sided")
#' str(stats)
#' 
#' # several unrelated contrasts
#' cls2 <- cls
#' cls[1:10] <- 1 # varying nx, ny
#' cls_mat <- cbind(cls, cls2)
#' stats <- rowWilcoxonTests(mat, categ = cls_mat, alternative = "two.sided")
#' str(stats)

rowWilcoxonTests <- function (mat, categ, alternative = c("two.sided", "less", "greater"), 
                              correct = TRUE) 
{
    alternative <- match.arg(alternative)
    stopifnot(all(categ %in% c(0, 1)))
    categ <- as.matrix(categ)
    levels(categ) <- NULL
    apply(categ, 2, categCheck, n = ncol(mat))
    B <- ncol(categ)
    m <- nrow(mat)
    n_obs <- nrow(categ)
    
    rks <- rowRanks(mat, ties.method = "average")
    nx <- colSums(categ)
    ny <- n_obs - nx
    
    min_stat <- nx * (nx + 1)/2 
    
    stats <- rks %*% categ
    # stats <- stats - min_stat 
    stats <- sweep(stats, MARGIN = 2, STATS = min_stat, FUN = "-")
    
    rks2 <- rks 
    mode(rks2) <- "integer"
    n_ties <- rowTabulates(rks2) 
    ties <- rowSums(n_ties^3 - n_ties) 
    # sigma <- sqrt((nx * ny/12) * ((ny + nx + 1) - ties/((ny + nx) * (ny + nx - 1))))
    # sigma <- t(sqrt((nx * ny/12) * ((ny + nx + 1) - matrix(rep(ties, B), nrow = B, byrow = T)/((ny + nx) * (ny + nx - 1)))))
    quotient <- sweep(matrix(rep(ties, B), ncol = B), MARGIN = 2, STATS = (ny + nx) * (ny + nx - 1), FUN = "/")
    difference <- sweep(-quotient, MARGIN = 2, STATS = (ny + nx + 1), FUN = "+")
    sigma <- sqrt(sweep(difference, MARGIN = 2, STATS = (nx * ny/12), FUN = "*"))
    
    # z <- stats - ny * nx/2
    z <- sweep(stats, MARGIN = 2, STATS = ny * nx/2, FUN = "-")
    CORRECTION <- 0
    if (correct) {
        CORRECTION <- switch(alternative, two.sided = sign(z) * 
                                 0.5, greater = 0.5, less = -0.5)
    }
    
    z <- (z - CORRECTION)/sigma
    p <- switch(alternative, 
                less = pnorm(z), 
                greater = pnorm(z, lower.tail = FALSE), 
                two.sided = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE)))
    
    est <- matrix(NA_real_, nrow = m, ncol = B)
    if (dim(categ)[2] == 1){
        wx <- which(categ == 1)
        est <- rowMedians(mat[, wx]) - rowMedians(mat[, -wx])
        stats <- as.vector(stats)
        p <- as.vector(p)
    } 
    list(p.value = p, statistic = stats, estimate = est)
}



# first version

rowWilcoxonTestsV1 <- function(mat, categ, 
                               alternative = c("two.sided", "less", "greater"), 
                               correct = TRUE) {
    alternative <- match.arg(alternative)
    stopifnot(all(categ %in% c(0, 1)))
    levels(categ) <- NULL
    
    if (is.vector(categ)) {
        return(rowWilcoxonTests1V1(mat, categ, 
                                   alternative = alternative, 
                                   correct = correct))
    } 
    stopifnot(is.matrix(categ))
    B <- ncol(categ)
    m <- nrow(mat)
    
    pval0 <- matrix(NA_real_, nrow = m, ncol = B) 
    stat0 <- matrix(NA_real_, nrow = m, ncol = B) 
    for (bb in 1:B) {
        rwt <- rowWilcoxonTests1V1(mat, categ[, bb], 
                                   alternative = alternative, 
                                   correct = correct)
        pval0[, bb] <- rwt$p.value
        stat0[, bb] <- rwt$statistic
    }
    
    list(p.value = pval0,
         statistic = stat0)
}

#' @importFrom matrixStats rowMedians
rowWilcoxonTests1V1 <- function(mat, categ, 
                                alternative = c("two.sided", "less", "greater"), 
                                correct = TRUE) {
    alternative <- match.arg(alternative)
    categCheck(categ, ncol(mat))
    
    rks <- rowRanks(mat, ties.method = "average")
    ny <- sum(categ == 0)
    nx <- sum(categ == 1)
    min_stat <- nx*(nx + 1)/2  ## mininmal rank sum
    
    wx <- which(categ == 1)
    stat <- rowSums(rks[, wx])
    stat <- stat - min_stat
    
    ## gaussian approximation for p-values (presence of ties or n > 50)
    ## source: 'stats:::wilcox.test.default'
    mode(rks) <- "integer"   # rowTabulates only takes integer values
    n_ties <- rowTabulates(rks)  
    ties <- rowSums(n_ties^3 - n_ties) 
    sigma <- sqrt((nx * ny / 12) * ((ny + nx + 1) - ties/((ny + nx) * (ny + nx - 1))))
    
    z <- stat - ny*nx/2
    CORRECTION <- 0
    if (correct) {
        CORRECTION <- switch(alternative, two.sided = sign(z) * 
                                 0.5, greater = 0.5, less = -0.5)
    }
    z <- (z - CORRECTION)/sigma
    p <- switch(alternative, 
                less = pnorm(z), 
                greater = pnorm(z, lower.tail = FALSE), 
                two.sided = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE)))
    
    est <- rowMedians(mat[, wx]) - rowMedians(mat[, -wx])
    list(statistic = stat,
         p.value = p, 
         estimate = est)
}


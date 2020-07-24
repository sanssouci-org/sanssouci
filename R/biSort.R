# bisort
# 
#  sorting a matrix decreasingly by row and then by columns
# 
#  The output matrix is sometimes referred to as "the Q matrix", following the
#  notation of Meinshausen (2006).
# 
#  @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of p-values 
#  under the null hypothesis. \describe{ \item{m}{is the number of
#  null hypotheses tested} \item{B}{is the number of Monte-Carlo samples}}
#  
#  @param kMax A scalar value between 1 and m that controls the depth of the
#  partial sorting of each row.
#  
#  @param Rcpp If \code{TRUE} (the default), sorting operations are performed
#  in C++.  Otherwise, the are performed using \code{R}'s \code{base::sort}.
#  
#  @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#  
#  @return A matrix of the same dimensions as the \code{mat} and
#  whose columns, then rows, have been sorted increasingly.  For
#  convenience the matrix where only columns are sorted is returned
#  as the attribute \code{'kmaxH0'} of \code{Q}.
#  
#  @details The output matrix is sometimes referred to as "the Q
#  matrix", following the notation of Meinshausen (2006).
#  @examples
# 
#  m <- 1023
#  B <- 1e3
#  pi0 <- 0.8
#  flavor <- c("independent", "equi-correlated", "3-factor model")[2]
#  rho <- 0.2
#  mat <- gaussianTestStatistics(m, B, dep = "equi", param = rho)$X0
#  pvals <- pnorm(mat, lower.tail = FALSE)
#  Q <- sansSouci:::bisort(mat, Rcpp=TRUE)
#  QR <- sansSouci:::bisort(mat, Rcpp=FALSE)
#  identical(Q,QR)
# 
#  kmaxH0 <- attr(Q, "kmaxH0")
# 
bisort <- function(mat, kMax=nrow(mat), Rcpp=TRUE) {
    m <- nrow(mat)
    B <- ncol(mat)
    stopifnot( 1 <= kMax && kMax <= m)
    ## k-max of the test statistics under H0:
    if (Rcpp) {
        kmaxH0 <- partialColSortDesc(mat, kMax);
        Q <- rowSortDesc(kmaxH0)
    } else {
        kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)
        kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]              ## truncate to [1, kMax]
        kmaxH0s <- apply(kmaxH0, 1, sort, decreasing=TRUE)
        Q <- t(kmaxH0s)
    }
    attr(Q, "kmaxH0") <- kmaxH0
    return(Q);
}


###########################################################################
## HISTORY:
#
# 2016-05-24
# o Created from getJoinFWERThresholds.
###########################################################################




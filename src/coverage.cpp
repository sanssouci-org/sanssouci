#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Empirical coverage of a threshold family
//
// @param thr A numeric vector of length \code{m}, a threshold family. Should
//   be sorted in decreasing order.
// @param Z A \code{m} x \code{B} matrix of \code{B} Monte-Carlo
//    samples of *ordered* test statistics under the null
//    hypothesis. Each column should be sorted in decreasing order.
// @return The empirical coverage of \code{thr} according to \code{Z}.
//
// @details The empirical coverage of \code{thr} according to
//   \code{Z} is defined as the empirical probability that there exists \eqn{k}
//   in \eqn{1 \dots m} such that \eqn{k}-th max of the joint distribution of
//   the test statistics under the null hypothesis exceeds the \eqn{k}-th
//   largest value of \code{thr}, that is, \code{thr[k]}.
//
// [[Rcpp::export]]
NumericVector empiricalCoverageO(NumericVector thr, arma::mat Z) {
    // Calculate P_n(k-max(thr) > thr(k))
    int m = thr.size();
    int B = Z.n_cols;
    int m2 = Z.n_rows;
    if (m != m2) {
        stop("Length of argument 'thr' should match number of rows of argument 'Z'");
    }
    NumericVector prob;
    int count=0;
    for (int bb=0; bb<B; bb++) {
        for (int ii=0; ii<m; ii++) {
            if (Z(ii,bb)> thr(ii)) {
                count = count + 1;
                break;
            }
        }
    }
    prob = count/(double)B;
    return prob;
}
/*** R
m <- 10
B <- 20
thr <- sort(rnorm(m), decreasing=TRUE)
Z <- matrix(rnorm(m*B), m, B)
p <- empiricalCoverageO(thr, Z)

## R version:
isAbove <- sweep(Z, 1, thr,  ">")
nAbove <- colSums(isAbove)
pR <- mean(nAbove>0)

stopifnot(identical(pR,p))
*/

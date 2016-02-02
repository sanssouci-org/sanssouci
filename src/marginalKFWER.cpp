#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector marginalKFWER(NumericVector thr, arma::mat Z) {
    // Calculate P_n(k-max(thr) > thr(k))
    int m = thr.size();
    int B = Z.n_cols;
    int m2 = Z.n_rows;
    if (m != m2) {
        stop("Length of argument 'thr' should match number of rows of argument 'Z'");
    }
    NumericVector probs(m);
    int count;
    for (int ii=0; ii<m; ii++) {
        count=0;
        double val = thr(ii);
        arma::rowvec z = Z.row(ii);
        for (int bb=0; bb<B; bb++) {
            if (z(bb)> val) {
             count = count+1;
            }
        }
        probs(ii) = count/(double)B;
    }
    return probs;
}

/*** R
m <- 10
B <- 20
thr <- sort(rnorm(m), decreasing=TRUE)
Z <- matrix(rnorm(m*B), m, B)
marginalKFWER(thr, Z)
*/

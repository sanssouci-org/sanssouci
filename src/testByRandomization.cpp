#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

arma::mat testByPermutation(arma::mat X, NumericVector cls, double B) {
    int m = X.n_rows;
    //    int n = X.n_cols;
    arma::mat T(m, B, arma::fill::zeros);
    
    // arma::vec eps;
    // arma::mat Tb;
    //    RNGScope scope;
    int bb, ii;
    for (bb=0; bb<B; bb++) {
        for (ii=0; ii<m; ii++) {
            // TODO!
        }
    }
    return(T);
}


/*** R
m <- 12
n <- 38
B <- 10
X <- matrix(rnorm(m*n), ncol=n)
cls <- rep(c(0, 1), times=c(27, n-27))

# set.seed(123)
# T <- testByPermutation(X, cls, B)
# 
# set.seed(123)
# TR <- testByPermutationR(X, cls, B)
# max(abs(T-TR))
*/

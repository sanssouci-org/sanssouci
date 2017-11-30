#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat testBySignFlipping(arma::mat X, double B) {
    int m = X.n_rows;
    int n = X.n_cols;
    arma::mat T(m, B, arma::fill::zeros);

    arma::vec eps;
    arma::mat Tb;
//    RNGScope scope;
    double isq = 1 / sqrt(n);
    int bb, ii;
    for (bb=0; bb<B; bb++) {
        eps = Rcpp::rbinom(n, 1, 0.5)*2 - 1;  // signs
        for (ii=0; ii<m; ii++) {
            Tb = X.row(ii) * eps;
            T(ii, bb) = Tb(0) * isq;
        }
    }
   return(T);
}


// [[Rcpp::export]]
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

set.seed(123)
T <- testBySignFlipping(X, B)

set.seed(123)
TR <- testBySignFlippingR(X, B)
max(abs(T-TR))

# set.seed(123)
# T <- testByPermutation(X, cls, B)
# 
# set.seed(123)
# TR <- testByPermutationR(X, cls, B)
# max(abs(T-TR))

*/

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat partialColSortDesc(arma::mat X, int k) {
// partial sorting by column in descending order
    int B = X.n_cols;
    arma::mat Y(k, B, arma::fill::zeros);
    for (int bb=0; bb<B; bb++) {
        arma::vec x = -X.col(bb);  // descending order
        std::nth_element(x.begin(), x.begin()+k, x.end());
        std::sort(x.begin(), x.begin()+k);
//        Y.col(bb) = -x;
        for (int rr=0; rr<k; rr++) {
            Y(rr,bb) = -x(rr);    // back to original values
        }
    }
    return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
B <- partialColSortDesc(A, 2);
print(B);

if (FALSE) {
    A <- matrix(rnorm(1e6), 10000, 1000);
    kmax <- 1000;
    benchmark(partialColSortDesc(A, kmax),
              -apply(-A, 2, sort, partial=kmax), replications=1)
}
*/

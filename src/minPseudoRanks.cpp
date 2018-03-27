#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Sorting of the rows of a matrix by descending order
//
// @param X A numeric matrix.
// @param Y A numeric matrix.
// @return A numeric matrix whose rows are sorted by descending order.
//
// @export
// [[Rcpp::export]]
Rcpp::NumericVector minPseudoRanks(arma::mat X, arma::mat Y) {
    int B = X.n_cols;
    int c = Y.n_rows;
    arma::vec z(B);
    arma::rowvec xk(B);
    double ykb;
    int ctr = 0;

    for (int bb=0; bb<B; bb++) {
      z(bb) = B; // initialize to the largest admissible value
    }
    for (int kk=0; kk<c; kk++) {
      xk = X.row(kk);
      for (int bb=0; bb<B; bb++) {
	ykb = Y(kk, bb);
	ctr = arma::sum(xk >= ykb);
	if (ctr < z(bb)) {
	  z(bb) = ctr;
	}
      }
    }
    z = z/B;

    // forget dimensions of 'x'
    Rcpp::NumericVector zz = Rcpp::wrap(z);
    zz.attr("dim") = R_NilValue;
    return zz;
}


/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
B <- rowSortDesc(A);
C <- rowSortDesc(A[1:2, ]);
D <- minPseudoRanks(B, C);
print(D);
*/

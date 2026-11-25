#' Mini batch version of a row_test_fun function.
#'
#' @param row_test_fun A function taking as inputs Y, categ and Alternative and
#' giving as a result a matrix of $p$-values.
#' Example \code{\link{row_welch_tests}} or \code{\link{row_wilcoxon_tests}}.
#' @param Y A \eqn{m \times n} numeric matrix whose rows correspond to variables
#' and columns to observations
#' @param categ Either a numeric vector of \code{n} categories in \eqn{0, 1} for
#'   the observations, or a \code{n x B} matrix stacking \code{B} such vectors
#'   (typically permutations of an original vector of size \code{n})
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less". As in
#'   \code{\link{t.test}}, alternative = "greater" is the alternative that class
#'   1 has a larger mean than class 0.
#' @param max_batch_size A number specifying the maximum size of information in
#'   each batch
#' @param m the number of hypotheses
#'
#' @return the p-values for the test
#' @export
#'
#' @examples
#' p <- 53
#' n <- 35
#' mat <- matrix(rnorm(p * n), ncol = n)
#' cls <- rep(c(0, 1), times = c(10, n - 10))
#' minibatch_row_test_fun(
#'   row_test_fun = row_welch_tests, Y = mat,
#'   categ = cls, alternative = "two.sided",
#'   max_batch_size = 1e5
#' )
minibatch_row_test_fun <- function(row_test_fun, Y, categ,
                                  alternative = c("two.sided", "less",
                                                  "greater"),
                                  max_batch_size = 1e6, m = nrow(Y)) {
  alternative <- match.arg(alternative)
  categ <- as.matrix(categ)

  B <- dim(categ)[2]
  nb_batch <- ceiling(B * m / max_batch_size)

  idxs <- rep(1:nb_batch, each = ceiling(B / nb_batch))[1:B]

  p_values <- matrix(NA_real_, nrow = m, ncol = B)
  for (batch_id in unique(idxs)) {
    id_batch <- which(idxs == batch_id)
    categ_batch <- categ[, id_batch]
    p_values[, id_batch] <- row_test_fun(Y,
      categ_batch,
      alternative = alternative
    )$p.value
  }
  return(p_values)
}

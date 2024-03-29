% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rowWilcoxonTests.R
\name{rowWilcoxonTests}
\alias{rowWilcoxonTests}
\title{Wilcoxon rank sum tests for each row of a matrix}
\usage{
rowWilcoxonTests(
  mat,
  categ,
  alternative = c("two.sided", "less", "greater"),
  correct = TRUE
)
}
\arguments{
\item{mat}{A \code{m x n} numeric matrix whose rows correspond to variables and columns to
observations}

\item{categ}{Either a numeric vector of \code{n} categories in \eqn{0, 1} for
the observations, or a \code{n x B} matrix stacking \code{B} such vectors
(typically permutations of an original vector of size \code{n})}

\item{alternative}{A character string specifying the alternative hypothesis.
Must be one of "two.sided" (default), "greater" or "less". As in
\code{\link{wilcox.test}}, alternative = "greater" is the alternative that
class 1 is shifted to the right of class 0.}

\item{correct}{A logical indicating whether to apply continuity correction in
the normal approximation for the p-value.}
}
\value{
A list containing the following components:
\describe{ \item{statistic}{the value of the statistics} \item{p.value}{the
p-values for the tests}}

A list containing the following components:
\describe{
\item{statistic}{the value of the statistics}
\item{p.value}{the p-values for the tests}
\item{estimate}{the median difference between groups (only calculated if \code{B=1} for computational efficiency)}}
Each of these elements is a matrix of size \code{m x B}, coerced to a vector of length \code{m} if \code{B=1}
}
\description{
Wilcoxon rank sum tests for each row of a matrix
}
\details{
This function performs \code{m x B} Wilcoxon T tests on
\code{n} observations. It is vectorized along the rows of \code{mat}. This
makes the code much faster than using loops of 'apply' functions,
especially for high-dimensional problems (small n and large m) because the
overhead of the call to the 'wilcox.test' function is avoided. Note that it
is not vectorized along the columns of \code{categ} (if any), as a basic
'for' loop is used.

The p-values are computed using the normal approximation as
described in the \code{\link{wilcox.test}} function. The exact p-values
(which can be useful for small samples with no ties) are not implemented
(yet).

For simplicity, 'estimate' returns the difference between the group medians, which does not match the component 'estimate' output by \code{wilcox.test}
}
\examples{

p <- 200
n <- 50
mat <- matrix(rnorm(p*n), ncol = n)
cls <- rep(c(0, 1), each = n/2)

stats <- rowWilcoxonTests(mat, categ = cls, alternative = "two.sided")
str(stats)

# permutation of class labels
cls_perm <- replicate(11, sample(cls))
stats <- rowWilcoxonTests(mat, categ = cls_perm, alternative = "two.sided")
str(stats)

# several unrelated contrasts
cls2 <- cls
cls[1:10] <- 1 # varying nx, ny
cls_mat <- cbind(cls, cls2)
stats <- rowWilcoxonTests(mat, categ = cls_mat, alternative = "two.sided")
str(stats)
}
\seealso{
wilcox.test
}
\author{
Gilles Blanchard, Pierre Neuvial and Etienne Roquain
}

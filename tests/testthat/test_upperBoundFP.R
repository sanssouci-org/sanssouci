context("Upper bound on the number of false positives")

test_that("curveMaxFP flavors give identical results", {
              m <- 13
              it <- c(1:2, 6:10, 13, 18:22)
              x <- sort(runif(2*m), decreasing=TRUE)
              T <- x[it]
              s <- x[-it]

              expect_equal(curveMaxFP(T, s, flavor="BNR2014"),
                           curveMaxFP(T, s, flavor="BNR2016"),
                           curveMaxFP(T, s, flavor="Mein2006"))

              ntests <- 10
              for (tt in 1:ntests) {
                  m <- 1e4
                  it <- sort(sample(2*m, m, replace=FALSE))
                  x <- sort(runif(2*m), decreasing=TRUE)
                  T <- x[it]
                  s <- x[-it]

                  expect_equal(curveMaxFP(T, s, flavor="BNR2016"),
                               curveMaxFP(T, s, flavor="Mein2006"))
              }
          })

test_that("curveMaxFP flavors give identical results with kMax", {
              m <- 13
              it <- c(1:2, 6:8, 13, 18)
              kMax <- length(it)
              x <- sort(runif(m+kMax), decreasing=TRUE)
              T <- x[it]  ## of length 'm'
              s <- x[-it] ## of length 'kMax'

              expect_equal(curveMaxFP(T, s, flavor="BNR2014"),
                           curveMaxFP(T, s, flavor="BNR2016"),
                           curveMaxFP(T, s, flavor="Mein2006"))

              ntests <- 10
              for (tt in 1:ntests) {
                  m <- 1e4
                  for (kMax in c(103, 207, 1011)) {
                      it <- sort(sample(m+kMax, m, replace=FALSE))
                      x <- sort(runif(m+kMax), decreasing=TRUE)
                      T <- x[it]
                      s <- x[-it]

                      expect_equal(curveMaxFP(T, s, flavor="BNR2016"),
                                   curveMaxFP(T, s, flavor="Mein2006"))
                  }
              }
          })

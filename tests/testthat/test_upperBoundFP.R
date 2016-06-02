context("Upper bound on the number of false positives")

test_that("upperBoundFP flavors give identical results", {
              m <- 13
              it <- c(1:2, 6:10, 13, 18:22)
              x <- sort(runif(2*m), decreasing=TRUE)
              T <- x[it]
              s <- x[-it]

              expect_equal(upperBoundFP(T, s, flavor="BNR2014"),
                           upperBoundFP(T, s, flavor="BNR2016"),
                           upperBoundFP(T, s, flavor="Mein2006"))

              ntests <- 10
              for (tt in 1:ntests) {
                  m <- 1e4
                  it <- sort(sample(2*m, m, replace=FALSE))
                  x <- sort(runif(2*m), decreasing=TRUE)
                  T <- x[it]
                  s <- x[-it]

                  expect_equal(upperBoundFP(T, s, flavor="BNR2016"),
                               upperBoundFP(T, s, flavor="Mein2006"))
              }
          })

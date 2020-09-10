context("Upper bound on the number of false positives")

test_that("curveMaxFP flavors give identical results for kMax = m", {
              m <- 13
              it <- c(1:2, 6:10, 13, 18:22)
              x <- sort(runif(2*m))
              p <- x[it]  
              thr <- x[-it]

              expect_equal(curveMaxFP(p, thr, flavor="BNR2014"),
                           curveMaxFP(p, thr, flavor="BNR2016"),
                           curveMaxFP(p, thr, flavor="Mein2006"))

              ntests <- 10
              for (tt in 1:ntests) {
                  m <- 1e3
                  it <- sort(sample(2*m, m, replace=FALSE))
                  x <- sort(runif(2*m))
                  p <- x[it]  
                  thr <- x[-it]

                  expect_equal(curveMaxFP(p, thr, flavor="BNR2016"),
                               curveMaxFP(p, thr, flavor="Mein2006"))
              }
          })

test_that("curveMaxFP flavors give identical results with kMax", {
              m <- 13
              it <- c(1:2, 6:8, 13, 18)
              kMax <- length(it)
              x <- sort(runif(m+kMax))
              p <- x[it]    ## of length 'm'
              thr <- x[-it] ## of length 'kMax'

              expect_equal(curveMaxFP(p, thr, flavor="BNR2014"),
                           curveMaxFP(p, thr, flavor="BNR2016"),
                           curveMaxFP(p, thr, flavor="Mein2006"))

              ntests <- 10
              for (tt in 1:ntests) {
                  m <- 20
                  for (kMax in c(13, 27, 91)) {
                      it <- sort(sample(m+kMax, m, replace=FALSE))
                      x <- sort(runif(m+kMax))
                      p <- x[it]  
                      thr <- x[-it]
    
                      expect_equal(curveMaxFP(p, thr, flavor="BNR2014"),
                                   curveMaxFP(p, thr, flavor="BNR2016"),
                                   curveMaxFP(p, thr, flavor="Mein2006"))
                  }
              }
          })

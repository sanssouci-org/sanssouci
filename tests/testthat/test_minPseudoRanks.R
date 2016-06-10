context("Calculation of the pivotal statistic of BNR for kFWER threshold family")

test_that("minPseudoRanks gives identical results using R or C++ implementation", {
              m <- 1000
              B <- 101
              m0 <- 800


              ntests <- 10
              for (tt in 1:ntests) {
                  mat <- simulateGaussianNullsFromFactorModel(m=m, n=B, flavor="equi-correlated", rho=0)
                  C <- sort(sample(m0))
                  c <- length(m0)
                  kmaxH0 <- partialColSortDesc(mat, c);  ## no need to go further than c!
                  kmaxH0C <- partialColSortDesc(mat[C, ], c);


                  expect_equal(minPseudoRanks(kmaxH0, kmaxH0C),
                               minPseudoRanksR(kmaxH0, kmaxH0C))
              }
          })

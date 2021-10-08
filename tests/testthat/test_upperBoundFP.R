context("Upper bound on the number of false positives")

test_that("curveMaxFP flavors give identical results for kMax = s", {
    m <- 13
    it <- c(1:2, 6:10, 13, 18:22)
    x <- sort(runif(2*m))
    p <- x[it]  
    thr <- x[-it]
    
    ub <- curveMaxFP_old(p, thr, flavor = "BNR2016")
    expect_equal(curveMaxFP_old(p, thr, flavor = "BNR2014"), ub)
    expect_equal(curveMaxFP_old(p, thr, flavor = "Mein2006"), ub)
    expect_equal(curveMaxFP_ECN(p, thr), ub)
    expect_equal(curveMaxFP(p, thr), ub)
    
    m <- 1e3
    it <- sort(sample(2*m, m, replace=FALSE))
    x <- sort(runif(2*m))
    p <- x[it]  
    thr <- x[-it]
    
    ub <- curveMaxFP_old(p, thr, flavor = "BNR2016")
    expect_equal(curveMaxFP_old(p, thr, flavor = "Mein2006"), ub)
    expect_equal(curveMaxFP_ECN(p, thr), ub)
    expect_equal(curveMaxFP(p, thr), ub)
})

test_that("curveMaxFP flavors give identical results with kMax < s", {
    m <- 20
    for (kMax in c(13, 27, 91)) {
        it <- sort(sample(m + kMax, m, replace=FALSE))
        x <- sort(runif(m + kMax))
        p <- x[it]  
        thr <- x[-it]
        
        ub <- curveMaxFP_old(p, thr, flavor = "BNR2016")
        expect_equal(curveMaxFP_old(p, thr, flavor = "BNR2014"), ub)
        expect_equal(curveMaxFP_old(p, thr, flavor = "Mein2006"), ub)
        expect_equal(curveMaxFP_ECN(p, thr), ub)
        expect_equal(curveMaxFP(p, thr), ub)
    }
    
    p <- c(1, 3, 5, 8, 9, 11, 12, 15, 16, 17, 18, 19)/20
    thr <- c(2, 4, 6, 7, 10, 13, 14)/20
    
    ub <- curveMaxFP_old(p, thr)
    expect_equal(curveMaxFP_ECN(p, thr), ub)
    expect_equal(curveMaxFP(p, thr), ub)
})

test_that("curveMaxFP flavors give identical results with kMax > s", {
    m <- 20
    for (kMax in c(13, 27, 91)) {
        it <- sort(sample(m + kMax, m, replace=FALSE))
        x <- sort(runif(m + kMax))
        p <- x[-it]  
        thr <- x[it]
        
        ub <- curveMaxFP_old(p, thr, flavor = "BNR2016")
        expect_equal(curveMaxFP_old(p, thr, flavor = "BNR2014"), ub)
        expect_equal(curveMaxFP_old(p, thr, flavor = "Mein2006"), ub)
        expect_equal(curveMaxFP_ECN(p, thr), ub)
        expect_equal(curveMaxFP(p, thr), ub)
    }
})

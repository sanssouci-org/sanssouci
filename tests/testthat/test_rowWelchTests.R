context("Calculation of Welch test statistics and p-values")

test_that("suffWelchTest gives identical results to welch.test", {

    p <- 1e3
    n <- 2e3
    mat <- matrix(rnorm(p*n), ncol=n)
    cls <- rep(c(0, 1), times=c(n/2, n-n/2))
    fwt <- rowWelchTests(mat, categ=cls)

    # Ordinary Welch t-test (reasonably fast for 1000 tests)
    wt <- apply(mat, 1, FUN=function(x) {
        tt <- t.test(x[cls == 1], x[cls == 0])  ## test stat positive if "1 > 0"
        c(tt$p.value, tt$statistic, tt$parameter)
    })
    dwt <- as.data.frame(t(wt))
    names(dwt) <- c("p.value", "statistic", "parameter")
    expect_equal(fwt$p.value, dwt$p.value, tolerance = 1e-12)
    expect_equal(fwt$statistic, dwt$statistic, tolerance = 1e-12)
    expect_equal(fwt$parameter, dwt$parameter, tolerance = 1e-12)
})


test_that("rowWelchTest gives identical results to welch.test", {

    mat <-     
    replicate(10, {
        x <- rnorm(1234)
        y <- rnorm(2345)
        target <- t.test(x, y)
        swt <- suffWelchTest(mean(x), mean(y), sd(x), sd(y), length(x), length(y))
        
        expect_equivalent(swt$statistic, target$statistic)
        expect_equal(swt$p.value, target$p.value)
        expect_equivalent(swt$parameter, target$parameter)
    })
})

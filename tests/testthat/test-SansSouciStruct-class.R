test_that("Correctness of the constructor of SansSouciStruct", {
    s <- 100
    q <- 7
    m <- s*2^q

    obj <- SansSouciDyadic(m, s, flavor = "window.size", direction = "top-down")
    expect_s3_class(obj, "SansSouciStruct")

    dd <- dyadic.from.window.size(m, s, method = 2)
    obj1 <- SansSouciStruct(dd$C, dd$leaf_list)
    expect_s3_class(obj, "SansSouciStruct")

    expect_identical(obj1, obj)
    rm(obj1)
    
    expect_equal(names(obj), c("input", "parameters", "output"))
    expect_null(obj$parameters)
    expect_null(obj$output)
    
    input <- obj$input
    expect_equal(names(input), c("struct", "leaves", "m"))
    expect_identical(input$struct, dd$C)
    expect_identical(input$leaves, dd$leaf_list)
    expect_length(input$leaves, 2^q)
    expect_length(input$struct, q+1)

    expect_equal(nHyp(obj), m)
    expect_equal(obj$input$m, m)
})


test_that("'fit.SansSouciStruct' reproduces the results of 'Vstar'", {
    s <- 100
    q <- 7
    m <- s*2^q
    
    obj <- SansSouciDyadic(m, s, flavor = "window.size", direction = "top-down")
    
    K1 <- 8
    r <- 0.9
    m1 <-  r*K1*s
    barmu <- 3
    
    mu <- gen.mu.leaves(m = m, K1 = K1, d = r, grouped = TRUE, setting = "const", 
                        barmu = barmu, leaf_list = obj$input$leaves)
    obj$input$truth <- as.numeric(mu != 0)
    pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
    
    ord <- order(pvalues)
    idxs <- round(c(seq(from = 1, to = 2*m1, length = 30),
                    seq(from = 2*m1, to = m, length = 10)[-1]))
    
    alpha <- 0.1

    # Oracle 1
    H0 <- which(mu == 0)
    V <- cumsum(ord %in% H0)
    V <- V[idxs]
    
    # Oracle 2
    res <- fit(obj, alpha = alpha, p.values = pvalues, family = "Oracle")
    expect_s3_class(res, "SansSouci")
    Vmat <- predict(res, what = "FP", all = TRUE)
    VV <- Vmat$bound[idxs]
    expect_identical(VV, V)
    
    # Simes 1
    Vmat <- confCurveFromFam(p.values = pvalues, refFamily = "Simes", 
                          param = alpha, what = "FP")
    V <- Vmat$bound[idxs]
    
    # Simes 2
    res_Simes <- fit(obj, alpha = alpha, p.values = pvalues, family = "Simes")
    Vmat <- predict(res_Simes, what = "FP", all = TRUE)
    VV <- Vmat$bound[idxs]
    expect_identical(VV, V)

    idxs1 <- head(idxs, 10)

    # DKWM 1 (tree)
    res_DKWM <- fit(obj, alpha, p.values = pvalues, 
                    family = "DKWM")
    V <- sapply(idxs1, FUN = function(ii) {
        predict(res_DKWM, S = seq_len(ii), what = "FP")
    })

    # DKWM 2 (tree)
    struct <- obj$input$struct
    leaves <- obj$input$leaves
    ZL <- zetas.tree(struct, leaves, zeta.DKWM, pvalues, alpha = alpha)
    VV <- sapply(idxs1, FUN = function(ii) {
        V.star(ord[1:ii], 
               C = struct, 
               ZL = ZL, 
               leaf_list = leaves)
    })
    expect_identical(VV, V)
    
    # DKWM 1 (part)
    res_DKWM <- fit(obj, alpha, p.values = pvalues, 
                    family = "DKWM", flavor = "partition")
    V <- sapply(idxs1, FUN = function(ii) {
        predict(res_DKWM, S = seq_len(ii), what = "FP")
    })
    
    # DKWM 2 (part)
    struct <- obj$input$struct
    CC <- struct[length(struct)]
    leaves <- obj$input$leaves
    ZL <- zetas.tree(CC, leaves, zeta.DKWM, pvalues, alpha = alpha)

    VV <- sapply(idxs1, FUN = function(ii) {
        V.star(ord[1:ii], 
               C = CC, 
               ZL = ZL, 
               leaf_list = leaves)
    })
    expect_identical(VV, V)
})
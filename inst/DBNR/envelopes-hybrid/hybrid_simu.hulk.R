simu.hulk <- function(m, 
                      s = 10, 
                      K1 = floor(m/(s * 4)), 
                      d = 1, 
                      barmu = 4,
                      setting = "const",
                      methods = c("tree", "part", "Simes", "Oracle"),
                      rho = 0, 
                      grouped = FALSE, 
                      random.leaves = FALSE, 
                      alphas = 0.05,
                      verbose = FALSE) {
    methods <- match.arg(methods, several.ok = TRUE)
    dd <- dyadic.from.window.size(m, s, method = 2)
    leaf_list <- dd$leaf_list
    mu <- gen.mu.leaves(m, K1, d, grouped, setting, barmu, leaf_list)
    # m1 <- sum(mu > 0)
    m1 <-  d*K1*s ## covers the case where barmu==0
    # xmax <- min(2 * m1, m)
    xmax <- min(4/3 * m1, m)
    #idxs1 <- c(1:max(xmax, s))
    idxs1 <- round(seq(from = 1, to = xmax, length = 10))
    idxs2 <- round(seq(from = xmax + 1, to = m, length = 10))
    idxs <- c(idxs1, idxs2)
    
    pvalues <- gen.p.values(m, mu, rho)
    #hf <- cherry::hommelFast(pvalues)
    C <- dd$C
    Cs <- list(tree = C,
               part = C[length(C)])
    rm(C)
    
    orders <- list(p.value = order(pvalues), 
              mu = order(mu, decreasing = TRUE))
    configs <- expand.grid(order = names(orders), method = methods, alpha = alphas)
    resList <- list()
    for (cc in 1:nrow(configs)) {
        if (verbose) print(configs[cc, ])
        ord <- configs[cc, "order"]
        oo <- orders[[ord]]
        meth <- configs[cc, "method"]
        alpha <- configs[cc, "alpha"]
        #print(configs[cc, ])
        if (verbose) print(meth)
        if (meth == "Oracle") {
            H0 <- which(mu == 0)  ## formally not true when barmu is 0...
            V <- cumsum(oo %in% H0)
            V <- V[idxs]
        } else if (meth == "Simes") {
            thr <- alpha * 1:m/m
            FP_Simes <- sansSouci:::curveMaxFP(pvalues, thr)
            V <- FP_Simes[idxs]
            # V2 <- idxs - sapply(idxs, FUN = function(ii) {
            #     posthocBySimes(pvalues, oo[1:ii], alpha) ## slow!
            # })
            # stopifnot(identical(V, V2))
        } else if (meth %in% c("tree", "part")) {
            ZL <- zetas.tree(Cs[[meth]], leaf_list, zeta.DKWM, pvalues, alpha = alpha)
            V <- sapply(idxs, FUN = function(ii) {
                V.star(oo[1:ii], Cs[[meth]], ZL, leaf_list)
            })
        } 
        res <- data.frame(idxs, value = V, method = meth, order = ord, alpha = alpha)
        resList[[cc]] <- res
    }
    resList
}

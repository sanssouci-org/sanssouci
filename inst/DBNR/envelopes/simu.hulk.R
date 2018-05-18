simu.hulk <- function(m, 
                      s = 10, 
                      K1 = floor(m/(s * 4)), 
                      d = 1, 
                      barmu = 4,
                      setting = "const",
                      methods = c("tree", "part", "Simes"),
                      rho = 0, 
                      grouped = FALSE, 
                      random.leaves = FALSE, 
                      alpha = 0.05,
                      SimesWeights = c(0.1, 0.2, 0.5, 0.8, 0.9, 1)) {
    dd <- dyadic.from.window.size(m, s, method = 2)
    leaf_list <- dd$leaf_list
    mu <- gen.mu.leaves(m, K1, d, grouped, setting, barmu, leaf_list)
    # m1 <- sum(mu > 0)
    m1 <-  d*K1*s ## covers the case where barmu==0
    range <- min(2 * m1, m)
    range <- min(4/3 * m1, m)
    idxs1 <- c(1:max(range, ss))
    idxs2 <- round(seq(from = max(range) + 1, to = m, length = 10))
    idxs <- c(idxs1, idxs2)
    
    pvalues <- gen.p.values(m, mu, rho)
    C <- dd$C
    Cs <- list(tree = C,
               part = C[length(C)])
    rm(C)
    
    orders <- list(p.value = order(pvalues), 
              mu = order(mu, decreasing = TRUE))
    configs <- expand.grid(order = names(orders), method = methods, SimesWeight = SimesWeights)
    resList <- list()
    for (cc in 1:nrow(configs)) {
        ord <- configs[cc, "order"]
        oo <- orders[[ord]]
        meth <- configs[cc, "method"]
        sw <- configs[cc, "SimesWeight"]
        print(cc)
        print(configs[cc, ])
        if (meth == "Simes") {
            V <- idxs - sapply(idxs, FUN = function(ii) {
                posthocBySimes(pvalues, oo[1:ii], alpha*sw)
            })
        } else {
            ZL <- zetas.tree(Cs[[meth]], leaf_list, zeta.DKWM, pvalues, alpha = alpha*sw)
            V <- sapply(idxs, FUN = function(ii) {
                V.star(oo[1:ii], Cs[[meth]], ZL, leaf_list)
            })
        }
        res <- data.frame(idxs, value = V, method = meth, order = ord, SimesWeight = sw)
        resList[[cc]] <- res
    }
    resList
}

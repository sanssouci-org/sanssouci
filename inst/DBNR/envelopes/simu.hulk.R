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
                      alpha = 0.05) {
    dd <- dyadic.from.window.size(m, s, method = 2)
    leaf_list <- dd$leaf_list
    mu <- gen.mu.leaves(m, K1, d, grouped, setting, barmu, leaf_list)
    # m1 <- sum(mu > 0)
    m1 <-  d*K1*s ## covers the case where barmu==0
    range <- min(2 * m1, m)
    idxs1 <- c(1:max(range, ss))
    idxs2 <- round(seq(from = max(range) + 1, to = m, length = 10))
    idxs <- c(idxs1, idxs2)
    
    pvalues <- gen.p.values(m, mu, rho)
    C <- dd$C
    Cs <- list(tree = C, 
               part = C[length(C)])
    rm(C)
    
    ZLs <- list(tree = zetas.tree(Cs[["tree"]], leaf_list, zeta.DKWM, pvalues, alpha = alpha),
                part = zetas.tree(Cs[["part"]], leaf_list, zeta.DKWM, pvalues, alpha = alpha))
    orders <- list(p.value = order(pvalues), 
              mu = order(mu, decreasing = TRUE))
    configs <- expand.grid(order = names(orders), method = methods)
    resList <- list()
    times <- numeric(0L)
    for (cc in 1:nrow(configs)) {
        ord <- configs[cc, "order"]
        oo <- orders[[ord]]
        meth <- configs[cc, "method"]
        tt <- system.time({
            VSimes <- function(alpha) {
                idxs - sapply(idxs, FUN = function(ii) {
                    posthocBySimes(pvalues, oo[1:ii], alpha)
                })
            }
            VDBNR <- function(alpha, meth) {
                ZL <- ZLs[[meth]]
                C <- Cs[[meth]]
                sapply(idxs, FUN = function(ii) {
                    V.star(oo[1:ii], C, ZL, leaf_list)
                })
            }
            if (meth == "Simes") {
                V <- VSimes(alpha)
            } else if (meth == "hybrid-0.5") {
                V <- pmin(VSimes(alpha/2), VDBNR(alpha/2, "tree"))
            } else if (meth == "hybrid-0.9") {
                V <- pmin(VSimes(alpha*0.9), VDBNR(alpha*0.1, "tree"))
            } else if (meth == "hybrid-0.1") {
                V <- pmin(VSimes(alpha*0.1), VDBNR(alpha*0.9, "tree"))
            } else {
                V <- VDBNR(alpha/2, meth)
            }
        })
        res <- data.frame(idxs, value = V, method = meth, order = ord)
        resList[[cc]] <- res
        times <- c(times, tt[["user.self"]])
    }
    names(times) <- paste(configs$order, configs$method, sep="-")
    attr(resList, "times") <- times
    resList
}

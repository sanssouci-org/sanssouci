gen.mu.noleaves <- function(m, pi0, barmu) {
    mu <- numeric(m)
    m0 <- floor(m * pi0)
    m1 <- m - m0
    mu[sample(seq(1, m), m1)] <- barmu
}

# n=30 plot(1:n,gauss_bloc(1,n),type='l')

gauss_bloc <- function(barmu, len) {
    # this choice of sigma2 is so that the values at the edge of the vector 'out' are = 1/200
    sigma2 <- len^2/(8 * log(200))
    out <- exp(-(seq_len(len) - floor(len/2))^2/(2 * sigma2))
    return(barmu * out/mean(out))
}

gen.mu.leaves <- function(m, K1, d, grouped, setting, barmu, leaf_list) {
    mu <- numeric(m)
    K <- length(leaf_list)
    if (K1 > K) 
        stop("K1>K,\nwe don't have so many leaves")
    active_leaves <- numeric(K1)
    if (grouped) {
        active_leaves <- seq(1, K1) + sample(seq(0, K - K1), 1)
    } else {
        active_leaves <- sample(seq(1, K), K1)
    }
    for (i in active_leaves) {
        length_leaf <- length(leaf_list[[i]])
        m1loc <- floor(length_leaf * d)
        signal <- sample(length_leaf, m1loc)
        mu[leaf_list[[i]][signal]] <- switch(setting, const = {
            barmu
        }, gauss = {
            gauss_bloc(barmu, length_leaf)[signal]
        }, poisson = {
            rpois(m1loc, 999 * barmu/1000) + barmu/1000
        })
    }
    return(mu)
}
gen.p.values <- function(m, mu, rho) {
    Z <- rnorm(m + 1, 0, 1)
    Y <- sqrt(1 - rho) * Z[seq_len(m)] + sqrt(rho) * Z[m + 1]
    return(1 - pnorm(Y + mu))
}
plotting <- function(C, ZL, leaf_list, method, pvalues, mu, alpha) {
    m <- length(pvalues)
    o <- order(pvalues)
    omu <- order(mu, decreasing = TRUE)
    vecVstar <- numeric(m)
    for (i in 1:m) {
        vecVstar[i] <- switch(method, threshold = {
            V.star(o[1:i], C, ZL, leaf_list)
        }, add = {
            V.star(omu[1:i], C, ZL, leaf_list)
        })
    }
    plot(1:m, vecVstar, type = "l")
    vequo <- numeric(m)
    for (i in 1:m) {
        vequo[i] <- switch(method, threshold = {
            min(sapply(1:m, function(k) min(sum(pvalues[o[1:i]] > alpha * k/m) + k - 1, i)))
        }, add = {
            min(sapply(1:m, function(k) min(sum(pvalues[omu[1:i]] > alpha * k/m) + k - 1, i)))
        })
    }
    lines(1:m, vequo, col = 2)
}

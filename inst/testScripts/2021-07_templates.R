library("sansSouci")
library("ggplot2")

# shape of linear vs beta templates
n_hyp <- 5e4
x <- 1:n_hyp/n_hyp
dat <- list(
    data.frame(x = x, y = t_inv_beta(x, n_hyp/10, n_hyp), 
               template = "Beta", k = n_hyp/10),
    data.frame(x = x, y = t_inv_linear(x, n_hyp/10, n_hyp), 
               template = "Linear", k = n_hyp/10),
    data.frame(x = x, y = t_inv_beta(x, 5*n_hyp/10, n_hyp), 
               template = "Beta", k = 5*n_hyp/10),
    data.frame(x = x, y = t_inv_linear(x, 5*n_hyp/10, n_hyp), 
               template = "Linear", k = 5*n_hyp/10),
    data.frame(x = x, y = t_inv_beta(x, 8*n_hyp/10, n_hyp), 
               template = "Beta", k = 8*n_hyp/10),
    data.frame(x = x, y = t_inv_linear(x, 8*n_hyp/10, n_hyp), 
               template = "Linear", k = 8*n_hyp/10))
dat <- Reduce(rbind, dat)
dat$k <- as.factor(dat$k)

p <- ggplot(dat, aes(x = x, y = y)) +
    geom_line(aes(color = k, linetype = template)) +
    ylim(c(0, 1.2)) +
    labs(x = expression(p[(k:H_0)]), y = expression(t_k(p[(k:H_0)]))) + 
    ggtitle(paste(n_hyp, "hypotheses")) + 
    geom_vline(xintercept = c(1, 5, 8)/10, linetype="dotted")
filename <- sprintf("inverse-templates_n_hyp=%s.png", n_hyp)
ggsave(p, filename = filename, scale = 0.75)

# A few realizations under indep and equi-correlation
B <- 100
rho <- 0.0
n_hyp <- 50000

set.seed(20210729)
obj <- SansSouciSim(m = n_hyp, rho = rho, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
res <- fit(obj, B = B, alpha = 0.1)
null_p <- res$output$p0

frac_k <- c(1, 5, 8)/10
seq_b <- sample(B, 6)

filename <- sprintf("beta-template_m=%s_rho=%s.pdf", n_hyp, rho)
pdf(filename, width = 4, height = 4)
title_ <- ifelse(rho == 0, "Independence", paste("Equicorrelation", rho))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = expression(p[(k:H_0)]), 
     ylab = expression(t_k(p[(k:H_0)])),
     main = title_)
for (k in frac_k*n_hyp)
    curve(t_inv_beta(x, k, n_hyp), add = TRUE)
for (i in seq(along = seq_b)) {
    b <- seq_b[i]
    x <- null_p[, b]
    sx <- sort(x)
    for (k in frac_k*n_hyp) {
        points(sx[k], t_inv_beta(sx[k], k, n_hyp), 
               col = i, pch = i)
    }
}
dev.off()
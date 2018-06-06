library("ggplot2")
m <- 12800;
s <- 1000;
K1 <- 8;
d <- 1;
m1 <-  d*K1*s;
barmu <- 3;

#Rprof("/tmp/hulk.Rout")
res <- simu.hulk(m = m, s = s, K1 = K1, grouped = TRUE, setting = "gauss", d = d, barmu = barmu, alphas = c(0.01, 0.05), verbose = TRUE) 
Rprof(NULL)
#summaryRprof("/tmp/hulk.Rout")

dat <- Reduce(rbind, res)
dat$level <- as.factor(dat$alpha)

xymax <- 4/3*m1;
ggplot(dat, aes(idxs, value, colour = method)) + 
    geom_line() +
    xlim(1, xymax) + ylim(0, xymax) +
    facet_grid(level~order, labeller = label_both) +
    ylab("Upper bound on the number of false positives") +
    xlab("sorted hypotheses")

ggplot(dat, aes(idxs, value, colour = method, linetype = level)) + 
    geom_line() +
    facet_wrap("order", labeller = label_both) +
    xlim(1, xymax) + ylim(0, xymax) +
    ylab("Upper bound on the number of false positives") +
    xlab("sorted hypotheses")

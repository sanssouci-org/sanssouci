library("ggplot2")
m <- 10000;
s <- 40;
K1 <- 8;
d <- 1;
barmu <- 3;
m1 <-  d*K1*s;
xymax <- 2*m1;

res <- simu.hulk(m = m, s = s, K1 = K1, grouped = TRUE, d = d, barmu = barmu) 
dat <- Reduce(rbind, res)
ggplot(dat, aes(idxs, value, colour = method, linetype = order)) + 
    geom_line() + 
    xlim(1, xymax) + ylim(0, xymax)

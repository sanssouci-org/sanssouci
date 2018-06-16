m <- 12800/4;
s <- 100;
K1 <- 8;
d <- 1;
m1 <-  d*K1*s;
barmu <- 3;
setting <- c("const", "gauss", "rgauss")[2]

#Rprof("/tmp/hulk.Rout")
res <- simu.hulk(m = m, s = s, K1 = K1, grouped = TRUE, setting = setting, d = d, barmu = barmu, alphas = c(0.0001, 0.05), verbose = TRUE) 
#Rprof(NULL)
#summaryRprof("/tmp/hulk.Rout")

## Vbar
dat <- Reduce(rbind, res)
dat$level <- as.factor(dat$alpha)

library("ggplot2")
lvls <- c("Oracle", "part", "Simes", "tree", "hybrid")
cols <- RColorBrewer::brewer.pal(length(lvls), "Set1")
names(cols) <- lvls

xymax <- 4/3*m1;
ggplot(dat, aes(idxs, value, colour = method)) + 
    geom_line() +
#    xlim(1, xymax) + ylim(0, xymax) +
    facet_grid(level~order, labeller = label_both) +
    ylab("Upper bound on the number of false positives") +
    xlab("sorted hypotheses") +
    scale_colour_manual(values = cols)

## Sbar
dat$S <- dat$idxs - dat$value

xmax <- 4/3*m1;
ymax <- max(dat$S);
ggplot(dat, aes(idxs, S, colour = method)) + 
    geom_line() +
    xlim(1, xmax) + ylim(0, ymax) +
    facet_grid(level~order, labeller = label_both) +
    ylab("Lower bound on the number of true positives") +
    xlab("sorted hypotheses")

ggplot(dat, aes(idxs, value, colour = method, linetype = level)) + 
    geom_line() +
    facet_wrap("order", labeller = label_both) +
    xlim(1, xymax) + ylim(0, xymax) +
    ylab("Upper bound on the number of false positives") +
    xlab("sorted hypotheses")

ggplot(dat, aes(idxs, S, colour = method, linetype = level)) + 
    geom_line() +
    facet_wrap("order", labeller = label_both) +
    xlim(1, xmax) + ylim(0, ymax) +
    ylab("Lower bound on the number of true positives") +
    xlab("sorted hypotheses")


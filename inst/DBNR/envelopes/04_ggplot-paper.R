resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

figPath <- R.utils::Arguments$getReadablePath("fig/DBNR")


grouped <- TRUE
setting <- "const"

filename <- sprintf("conf-env-alpha_grouped=%s_setting=%s.rds", grouped, setting)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
rm(pathname)
dim(dat)
        
library("magrittr")
library("dplyr")

gdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5 & barmu <= 4) %>%
    rename(r = d) %>%
    group_by(idxs, method, r, barmu, alpha, order) %>% 
    summarise(value = mean(value))

library("ggplot2")
xymax <- 800

## 1- overview plot for a fixed alpha
plotConfigs <- expand.grid(alpha = unique(gdat$alpha),
                           order = unique(gdat$order))

for (pp in 1:nrow(plotConfigs)) {
    aa <- plotConfigs[pp, "alpha"]
    oo <- plotConfigs[pp, "order"]
    ftag <- sprintf("_alpha=%s_order=%s", aa, oo)
    ftag <- gsub("\\.", "-", ftag)
    ftag <- sprintf("%s.pdf", ftag)
    ff <- gsub("\\.rds$", ftag, filename)
    pathname <- file.path(figPath, ff)
    
    sdat <- subset(gdat, alpha == aa & order == oo)
    ggplot(sdat, aes(idxs, value, colour = method)) + 
        geom_line() + 
        facet_grid(r ~ barmu, labeller = label_both) +
        xlim(1, xymax) + ylim(0, xymax) +   
        ylab(sprintf("%s%% upper confidence bound on the number of false positives", (1-aa)*100)) +
        xlab(sprintf("Hypotheses sorted by %s", oo)) +
        geom_abline(slope = 1, intercept = 0, col = "lightgray", lty = 2)
    ggsave(filename = pathname)
}

## 2- influence of alpha
plotConfigs <- expand.grid(r = unique(gdat$r),
                           barmu = unique(gdat$barmu),
                           order = unique(gdat$order))
for (pp in 1:nrow(plotConfigs)) {
    rr <- plotConfigs[pp, "r"]
    mu <- plotConfigs[pp, "barmu"]
    oo <- plotConfigs[pp, "order"]
    
    ftag <- sprintf("_r=%s_mu=%s_order=%s", rr, mu, oo)
    ftag <- gsub("\\.", "-", ftag)
    ftag <- sprintf("%s.pdf", ftag)
    ff <- gsub("\\.rds$", ftag, filename)
    pathname <- file.path(figPath, ff)
    
    sdat <- subset(gdat, r == rr & barmu == mu & order == oo & method != "Oracle")
    sdat$alpha <- as.factor(sdat$alpha)
    ggplot(sdat, aes(idxs, value, colour = method, linetype = alpha)) + 
        geom_line() + facet_grid(order ~ method, labeller = label_both) +
        xlim(1, xymax) + ylim(0, xymax) +
        ylab("Upper confidence bounds on the number of false positives") +
        xlab(sprintf("Hypotheses sorted by %s", oo))
    ggsave(filename = pathname, height = 4)
}

## 2b- influence of alpha -- all in one
plotConfigs <- expand.grid(order = unique(gdat$order))
for (pp in 1:nrow(plotConfigs)) {
    oo <- plotConfigs[pp, "order"]
    
    ftag <- sprintf("_order=%s", oo)
    ftag <- gsub("\\.", "-", ftag)
    ftag <- sprintf("%s.pdf", ftag)
    ff <- gsub("\\.rds$", ftag, filename)
    pathname <- file.path(figPath, ff)
    
    sdat <- subset(gdat, order == oo & barmu > 2 & r < 1)
    sdat$alpha <- as.factor(sdat$alpha)
    ggplot(sdat, aes(idxs, value, colour = method, linetype = alpha)) + 
        geom_line() + facet_grid(r ~ barmu, labeller = label_both) +
        xlim(1, xymax) + ylim(0, xymax) +
        ylab("Upper confidence bounds on the number of false positives") +
        xlab(sprintf("Hypotheses sorted by %s", oo)) +
        geom_abline(slope = 1, intercept = 0, col = "lightgray", lty = 2)
    ggsave(filename = pathname)
}
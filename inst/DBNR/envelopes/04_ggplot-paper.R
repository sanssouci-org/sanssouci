resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

figPath <- R.utils::Arguments$getReadablePath("fig/DBNR")

# line colors
lvls <- c("Oracle", "part", "Simes", "tree", "hybrid")
cols <- RColorBrewer::brewer.pal(length(lvls), "Set1")
names(cols) <- lvls
lblr <- label_bquote(
    rows = r == .(r), 
    cols = bar(mu) == .(barmu))

#filename <- sprintf("conf-env-alpha-hybrid_grouped=%s_setting=%s.rds", grouped, setting)
#filename <- sprintf("conf-env-alpha_grouped=%s_setting=%s.rds", grouped, setting)
figName <- "all-conf-env"
filename <- sprintf("2018-06-15_%s.rds", figName)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
rm(pathname)
dim(dat)
        
library("magrittr")
library("dplyr")
library("ggplot2")

xymax <- 800
ymax <- 601

## 1- overview plot for a fixed alpha
gdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5 & barmu <= 4 & grouped == TRUE 
               & alpha %in% c(0.0001, 0.001, 0.05)) %>%
    rename(r = d) %>%
    group_by(idxs, method, r, barmu, alpha, order, setting, grouped) %>% 
    summarise(value = mean(value))

plotConfigs <- expand.grid(alpha = unique(gdat$alpha),
                           order = unique(gdat$order),
                           setting = unique(gdat$setting))

for (pp in 1:nrow(plotConfigs)) {
    aa <- plotConfigs[pp, "alpha"]
    oo <- plotConfigs[pp, "order"]
    ss <- plotConfigs[pp, "setting"]
    ftag <- sprintf("alpha=%s_order=%s_setting=%s", aa, oo, ss)
    ftag <- gsub("\\.", "-", ftag)
    print(ftag)
    ff <- sprintf("%s_%s.pdf", figName, ftag)
    pathname <- file.path(figPath, ff)
    
    sdat <- subset(gdat, alpha == aa & order == oo & setting == ss)
    sdat$method <- factor(sdat$method, levels = lvls)
    
    ggplot(sdat, aes(idxs, idxs - value, colour = method)) + 
        geom_line() + 
        facet_grid(r ~ barmu, labeller = lblr) +
        # xlim(1, xymax) + ylim(0, xymax) +   
        # ylab(sprintf("%s%% upper confidence envelope on the number of false positives", (1-aa)*100)) +
        xlim(1, xymax) + ylim(0, ymax) +   
        ylab(sprintf("%s%% lower confidence envelope on the number of true positives", (1-aa)*100)) +
        xlab(sprintf("Hypotheses sorted by %s", oo)) +
        # geom_abline(slope = 1, intercept = 0, col = "lightgray", lty = 2) +
        scale_linetype(name = expression(alpha)) +
        scale_colour_manual(values = cols)
    ggsave(filename = pathname)
}


## 2- influence of alpha -- all in one
## now focus on order = p and grouped
gdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5 & barmu <= 4 & grouped == TRUE) %>%
    rename(r = d) %>%
    group_by(idxs, method, r, barmu, alpha, order, setting, grouped) %>% 
    summarise(value = mean(value))

plotConfigs <- expand.grid(
    setting = unique(gdat$setting),
    order = unique(gdat$order)
)

for (pp in 1:nrow(plotConfigs)) {
    oo <- plotConfigs[pp, "order"]
    ss <- plotConfigs[pp, "setting"]
    
    ftag <- sprintf("all-alpha_order=%s_setting=%s", oo, ss)
    ftag <- gsub("\\.", "-", ftag)
    ff <- sprintf("%s_%s.pdf", figName, ftag)
    pathname <- file.path(figPath, ff)
    
    sdat <- subset(gdat, setting == ss & order == oo & barmu > 2 & r < 1)
    sdat$alpha <- as.factor(sdat$alpha)
    ggplot(sdat, aes(idxs, idxs - value, colour = method, linetype = alpha)) + 
        geom_line() + facet_grid(r ~ barmu, labeller = lblr) +
        # xlim(1, xymax) + ylim(0, xymax) +
        # ylab("Upper confidence envelope on the number of false positives") +
        xlim(1, xymax) + ylim(0, ymax) +   
        ylab("Lower confidence envelope on the number of true positives") +
        xlab(sprintf("Hypotheses sorted by %s", oo)) +
        # geom_abline(slope = 1, intercept = 0, col = "lightgray", lty = 2) +
        scale_linetype(name = expression(alpha)) +
        scale_colour_manual(values = cols)
    ggsave(filename = pathname)
}

## 3- hybrid procedures!

gdat <- filter(dat, s == 100 & K1 == 8 
               & d > 0.5  & d < 1
               & barmu <= 4 & barmu > 2) %>%
    filter(grouped == TRUE) %>%
    filter(alpha %in% c(0.001, 0.049, 0.05)) %>%
    rename(r = d) %>%
    arrange(idxs, method, r, barmu, alpha, order)
simes <- filter(gdat, method == "Simes" & alpha == 0.049) 
tree <- filter(gdat, method == "tree" & alpha == 0.001)
bound <- pmin(select(simes, value), select(tree, value))
hyb <- simes %>% select(-value, -alpha, -method) %>%
    mutate(method = "hybrid", 
           alpha = 0.05, 
           value = bound$value)
gdat <- gdat %>% mutate(method = as.character(method))%>%
    bind_rows(hyb) %>%
    group_by(idxs, method, r, barmu, alpha, order, setting) %>% 
    summarise(value = mean(value))

plotConfigs <- expand.grid(
    setting = unique(ggdat$setting),
    order = unique(gdat$order)
)

for (pp in 1:nrow(plotConfigs)) {
    oo <- plotConfigs[pp, "order"]
    ss <- plotConfigs[pp, "setting"]
    
    ftag <- sprintf("hybrid_order=%s_setting=%s", oo, ss)
    ftag <- gsub("\\.", "-", ftag)
    ff <- sprintf("%s_%s.pdf", figName, ftag)
    pathname <- file.path(figPath, ff)

    sdat <- subset(gdat, order == oo & barmu >= 2 & r <= 1 & method != "part" & setting == ss)
    ssdat <- subset(sdat,  method == "Simes" & alpha %in% c(0.001, 0.05))
    tsdat <- subset(sdat,  method == "tree" & alpha %in% c(0.001, 0.05))
    hsdat <- subset(sdat,  method == "hybrid" & alpha == 0.05)
    osdat <- subset(sdat,  method == "Oracle" & alpha == 0.05)
    sdat <- rbind(ssdat, tsdat, hsdat, osdat)
    sdat$alpha <- as.factor(sdat$alpha)
    sdat$alpha <- factor(sdat$alpha, levels = rev(levels(sdat$alpha)))
    sdat$method <- factor(sdat$method, levels = lvls)
    
    ggplot(sdat, aes(idxs, idxs - value, colour = method, linetype = alpha)) + 
        geom_line() + 
        facet_grid(r ~ barmu, 
                   labeller = label_bquote(
                       rows = r == .(r), 
                       cols = bar(mu) == .(barmu))) +
        # xlim(1, xymax) + ylim(0, xymax) +   
        # ylab("Upper confidence envelope on the number of false positives") +
        xlim(1, xymax) + ylim(0, ymax) +
        ylab("Lower confidence envelope on the number of true positives") +
        xlab(sprintf("Hypotheses sorted by %s", oo)) +
        geom_abline(slope = 1, intercept = 0, col = "lightgray", lty = 2) +
        scale_linetype(name=expression(alpha)) +
        scale_colour_manual(values = cols)
    ggsave(filename = pathname)
}

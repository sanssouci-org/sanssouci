resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

figPath <- R.utils::Arguments$getReadablePath("fig/DBNR")


configs <- expand.grid(
    grouped = groupeds,
    setting = settings,
    stringsAsFactors = FALSE)

grouped <- TRUE
setting <- "const"

for (grouped in unique(configs$grouped)) {
    for (setting in unique (configs$setting)) {

        filename <- sprintf("conf-env-alpha_grouped=%s_setting=%s.rds", grouped, setting)
        pathname <- file.path(resPath, filename)
        dat <- readRDS(pathname)
        rm(pathname)
        dim(dat)
        
        names(dat)
        oo <- "p.value"
        sdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5 & order == oo)
        
        dim(sdat)
        
        library("magrittr")
        library("dplyr")
        gdat <- sdat %>% group_by(idxs, method, d, barmu, alpha) 
        pdat <- gdat %>% summarise(value = mean(value))
        
        library("ggplot2")
        m1 <-  unique(sdat$s)*unique(sdat$K1)*max(sdat$d)
        xymax <- 4/3*m1
        
        ff <- gsub("\\.rds$", ".pdf", filename)
        pathname <- file.path(figPath, ff)

        pdat <- subset(pdat, alpha %in% c(0.001, 0.01, 0.1, 0.25, 0.5))
        pdat$alpha <- as.factor(pdat$alpha)
        ggplot(pdat, aes(idxs, value, colour = method, linetype = alpha)) + 
            geom_line() + facet_grid(d ~ barmu, labeller = label_both) +
            xlim(1, xymax) + ylim(0, xymax) +
            ylab("Upper bound on the number of false positives") +
            xlab(sprintf("Hypotheses sorted by %s", oo))
        ggsave(filename = pathname)
    }
}
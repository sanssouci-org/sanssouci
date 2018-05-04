resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

figPath <- R.utils::Arguments$getReadablePath("fig/DBNR")


configs <- expand.grid(
    grouped = groupeds,
    setting = settings,
    stringsAsFactors = FALSE)

# grouped <- TRUE
# setting <- "const"

for (grouped in unique(configs$grouped)) {
    for (setting in unique (configs$setting)) {

        filename <- sprintf("conf-env_grouped=%s_setting=%s.rds", grouped, setting)
        pathname <- file.path(resPath, filename)
        dat <- readRDS(pathname)
        rm(pathname)
        dim(dat)
        
        names(dat)
        
        sdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5)
        dim(sdat)
        
        library("magrittr")
        library("dplyr")
        gdat <- sdat %>% group_by(idxs, method, order, d, barmu) 
        pdat <- gdat %>% summarise(value = mean(value))
        
        library("ggplot2")
        m1 <-  unique(sdat$s)*unique(sdat$K1)*max(sdat$d)
        xymax <- 2*m1
        
        ff <- gsub("\\.rds$", ".pdf", filename)
        pathname <- file.path(figPath, ff)
        
        ggplot(pdat, aes(idxs, value, colour = method, linetype = order)) + 
            geom_line() + facet_grid(d ~ barmu) +
            xlim(1, xymax) + ylim(0, xymax)
        ggsave(filename = pathname)
        
    }
}
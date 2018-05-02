resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

# configs <- expand.grid(
#     grouped = groupeds,
#     setting = settings, 
#     stringsAsFactors = FALSE)
# 
# for (grouped in unique(configs$grouped)) {
#     for (setting in unique (configs$setting)) {

grouped <- TRUE
setting <- "const"

filename <- sprintf("conf-env_grouped=%s_setting=%s.rds", grouped, setting)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
dim(dat)



library("ggplot2")
names(dat)


sdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5)

m1 <-  unique(sdat$s)*unique(sdat$K1)*max(sdat$d)
xymax <- 2*m1
                   
ggplot(sdat, aes(idxs, value, colour = method, linetype = order)) + 
    geom_line() + facet_grid(d ~ barmu) +
    xlim(1, xymax) + ylim(0, xymax)

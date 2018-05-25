resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getReadablePath(resPath)

resPath2 <- file.path(resPath, Sys.Date())
resPath2 <- R.utils::Arguments$getWritablePath(resPath2)

confs <- subset(configs, grouped & setting == "const")

## ms <- unique(configs$m)
## ss <- unique(configs$s)
## stopifnot(length(ms)==1)
## stopifnot(length(ss)==1)

for (grp in unique(configs$grouped)) {
    message("grouped=", grp)
    for (st in unique (configs$setting)) {
        message("setting=", st)
        confs <- subset(configs, grouped==grp & setting == st)
        str(confs)
        datList <- NULL
        for (cc in 1:nrow(confs)) {
            conf <- confs[cc, ]
            stag <- paste(
                "m=", conf[["m"]], "_",
                "s=", conf[["s"]], "_",
#                "m=", ms, "_",
#                "s=", ss, "_",
                "K1=", conf[["K1"]], "_",
                "d=", conf[["d"]], "_",
                "barmu=", conf[["barmu"]], "_",
                "grouped=", conf[["grouped"]], "_",
                "setting=", conf[["setting"]], sep = "")
            filename <- sprintf("conf-env_%s.rds", stag)
            pathname <- file.path(resPath, filename)
            datList[[cc]] <- readRDS(pathname)
        }
        dat <- Reduce(rbind, datList)
        print(nrow(dat))
        filename <- sprintf("conf-env_grouped=%s_setting=%s.rds", grp, st)
        pathname <- file.path(resPath2, filename)
        saveRDS(dat, pathname)
        rm(pathname)
    }
}

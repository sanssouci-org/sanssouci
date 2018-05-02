resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getWritablePath(resPath)

confs <- subset(configs, grouped & setting == "const")

for (grouped in unique(configs$grouped)) {
    for (setting in unique (configs$setting)) {
        print(setting)
        confs <- subset(configs, grouped & setting == "const")
        datList <- NULL
        for (cc in 1:nrow(confs)) {
            conf <- confs[cc, ]
            stag <- paste("m=", conf[["m"]], "_",
                          "s=", conf[["s"]], "_",
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
        filename <- sprintf("conf-env_grouped=%s_setting=%s.rds", grouped, setting)
        pathname <- file.path(resPath, filename)
        saveRDS(dat, pathname)
        rm(pathname)
    }
}

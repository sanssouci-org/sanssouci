library("future")
plan(multiprocess, workers = 3)
#plan(multiprocess, workers = 100)

resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getWritablePath(resPath)

nb <- nrow(configs)

for (cc in 1:nb) {
    conf <- configs[cc, ]
    stag <- paste("m=", conf[["m"]], "_",
                  "s=", conf[["s"]], "_",
                  "K1=", conf[["K1"]], "_",
                  "d=", conf[["d"]], "_",
                  "barmu=", conf[["barmu"]], "_",
                  "grouped=", conf[["grouped"]], "_",
                  "setting=", conf[["setting"]], sep = "")
    filename <- sprintf("conf_env-%s.rds", stag)
    dummy %<-% {
        res <- simu.hulk(m = conf[["m"]], 
                         s = conf[["s"]], 
                         K1 = conf[["K1"]], 
                         d = conf[["d"]], 
                         barmu = conf[["barmu"]],
                         grouped = conf[["grouped"]], 
                         setting = conf[["setting"]])
        dat <- Reduce(rbind, res)
        rownames(conf) <- NULL
        dat <- cbind(dat, conf)
        pathname <- file.path(resPath, filename)
        saveRDS(dat, pathname)
        ## dat <- readRDS(pathname)
    }
}

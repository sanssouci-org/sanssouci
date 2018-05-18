library("future")
plan(multiprocess, workers = 100)
#plan(multiprocess, workers = 100)

resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getWritablePath(resPath)

nb <- nrow(configs)

for (cc in 1:nb) {
    print(cc)
    conf <- configs[cc, ]
    stag <- paste("m=", conf[["m"]], "_",
                  "s=", conf[["s"]], "_",
                  "K1=", conf[["K1"]], "_",
                  "d=", conf[["d"]], "_",
                  "barmu=", conf[["barmu"]], "_",
                  "grouped=", conf[["grouped"]], "_",
                  "setting=", conf[["setting"]], sep = "")
    filename <- sprintf("conf-env_%s.rds", stag)
    dummy %<-% {
        sdatList <- list()
        for (rr in 1:repl) {
            res <- simu.hulk(m = conf[["m"]], 
                             s = conf[["s"]], 
                             K1 = conf[["K1"]], 
                             d = conf[["d"]], 
                             barmu = conf[["barmu"]],
                             grouped = conf[["grouped"]], 
                             setting = conf[["setting"]])
            sdat <- Reduce(rbind, res)
            sdat$replication <- rr
            sdatList[[rr]] <- sdat
        }
        dat <- Reduce(rbind, sdatList)
        
        rownames(conf) <- NULL
        dat <- cbind(dat, conf)
        pathname <- file.path(resPath, filename)
        saveRDS(dat, pathname)
        ## dat <- readRDS(pathname)
    }
}

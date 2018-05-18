library("future")
plan(multiprocess, workers = 100)
#plan(multiprocess, workers = 100)

resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getWritablePath(resPath)

nb <- nrow(configs)
#nb <- 40
#cc <- 220 ## (a nice one)
# > conf
# m   s K1   d barmu grouped setting
# 220 10000 100  8 0.9     4    TRUE   const
for (cc in 1:nb) {
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
        res <- simu.hulk(m = conf[["m"]], 
                         s = conf[["s"]], 
                         K1 = conf[["K1"]], 
                         d = conf[["d"]], 
                         barmu = conf[["barmu"]],
                         grouped = conf[["grouped"]], 
                         setting = conf[["setting"]],
                         methods = c("tree", "part", "Simes", "hybrid-0.5", "hybrid-0.9", "hybrid-0.1"))
        dat <- Reduce(rbind, res)
        rownames(conf) <- NULL
        dat <- cbind(dat, conf)
        pathname <- file.path(resPath, filename)
        saveRDS(dat, pathname)
        ## dat <- readRDS(pathname)
    }
}

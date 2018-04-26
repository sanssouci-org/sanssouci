resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- R.utils::Arguments$getWritablePath(resPath)

ms <- 10000
ss <- c(10, 20, 50, 100)
ds <- c(0, 0.5, 0.75, 0.9, 1)
barmus <- c(2, 3, 4, 5)
K1s <- c(1, 4, 8, 16)
groupeds <- c(TRUE, FALSE)
settings <- c("const", "gauss", "poisson")

configs <- expand.grid(
    m = ms,
    s = ss,
    K1 = K1s,
    d = ds,
    barmu = barmus,
    grouped = groupeds,
    setting = settings, 
    stringsAsFactors = FALSE)

nb <- nrow(configs)
nb <- 2

for (cc in 1:nb) {
    conf <- configs[cc, ]
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
    stag <- paste("m=", conf[["m"]], "_",
                  "s=", conf[["s"]], "_",
                  "K1=", conf[["K1"]], "_",
                  "d=", conf[["d"]], "_",
                  "barmu=", conf[["barmu"]], "_",
                  "grouped=", conf[["grouped"]], "_",
                  "setting=", conf[["setting"]], sep = "")
    filename <- sprintf("conf_env-%s.rds", stag)
    pathname <- file.path(resPath, filename)
    saveRDS(dat, pathname)
}



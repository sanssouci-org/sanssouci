resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- file.path(resPath, "2018-06-06")
resPath <- R.utils::Arguments$getWritablePath(resPath)

## confs <- subset(configs, grouped & setting == "const")

## ms <- unique(configs$m)
## ss <- unique(configs$s)
## stopifnot(length(ms)==1)
## stopifnot(length(ss)==1)

## ## one big file!
##for (grp in unique(configs$grouped)) {
##    message("grouped=", grp)
##    for (st in unique (configs$setting)) {
##        message("setting=", st)
##        confs <- subset(configs, grouped==grp & setting == st)
filename <- sprintf("all-conf-env-rgauss.txt")
opathname <- file.path(resPath, filename)
confs <- configs
str(confs)

for (cc in 1:nrow(confs)) {
#for (cc in 1:3) {
    conf <- confs[cc, ]
    if (cc %% 100 == 1) {
        print(conf)
    }
    stag <- paste(
        "m=", conf[["m"]], "_",
        "s=", conf[["s"]], "_",
        "K1=", conf[["K1"]], "_",
        "d=", conf[["d"]], "_",
        "barmu=", conf[["barmu"]], "_",
        "grouped=", conf[["grouped"]], "_",
        "setting=", conf[["setting"]], sep = "")
    filename <- sprintf("conf-env_%s.rds", stag)
    pathname <- file.path(resPath, filename)
    if (!file.exists(pathname)) {
        print(pathname)
    } else {
        dat <- readRDS(pathname)
        write.table(dat, file = opathname, quote = FALSE,
                    append = (cc>1), row.names = FALSE,
                    col.names = (cc==1))
    }
}
dat <- read.table(opathname, header = TRUE)
print(nrow(dat))
filename <- sprintf("all-conf-env-rgauss.rds")
pathname <- file.path(resPath, filename)
saveRDS(dat, pathname)
rm(pathname, opathname)


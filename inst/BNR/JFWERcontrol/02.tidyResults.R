## source("inst/BNR/averagePower/00.setup.R")

rpath <- file.path("~/Documents/Packages/sanssouci", path)

molten <- FALSE

if (molten) {
    fls <- list.files(rpath, pattern=",molten", full.names=TRUE)
    id <- gsub(",molten.rds$", "", basename(fls))
    names(fls) <- id
    dat <- plyr::ldply(fls, readRDS, .id="id")
    names(dat)[match("dep", names(dat))] <- "rho"
} else {
    fls <- list.files(rpath, full.names=TRUE)
    fls <- fls[-grep("molten", fls)]
    id <- gsub(".rds$", "", basename(fls))
    names(fls) <- id
    dat <- plyr::ldply(fls, readRDS, .id="id")
    names(dat)[match("dep", names(dat))] <- "rho"
}

head(dat)
filename <- sprintf("%s,%s.rds", sname, pname)
pathname <- file.path(resPath, filename)
#saveRDS(dat, file=pathname)

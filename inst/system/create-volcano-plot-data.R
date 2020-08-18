dataSets <- c("chiaretti", "golub", "hedenfalk", "loi", "rosenwald")
dataSets <- dataSets[-4] ## Dropping Loi: too weak signal

getStats <- function(dataSet) {
    filename <- sprintf("%s.rds", dataSet)
    pathname <- file.path("local/data", filename)
    dat <- readRDS(pathname)
    rm(pathname)
    
    cn <- colnames(dat)[1]
    categ <- ifelse(colnames(dat) == cn, 1, 0) # map to 0/1
    res <- rowWelchTests(dat, categ)
    mres <- Reduce(cbind, res)
    colnames(mres) <- names(res)
    
    m <- nrow(mres)
    gid <- rownames(mres)
    if (is.null(gid)) {
        patt <- sprintf("gene_%%0%ds", floor(log10(m)) + 1)
        gid <- sprintf(patt, 1:m)
    }
    p <- mres[ ,"p.value"]
    adjp <- p.adjust(p, method = "BH")
    dat <- data.frame(id = gid, 
                 dataSet = dataSet, 
                 meanDiff = round(mres[, "meanDiff"], 3), 
                 p.value = p,
                 adjp.BH = round(adjp, 3),
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    dat
}

dfList <- lapply(dataSets, getStats)
volcano <- Reduce(rbind, dfList)
save(volcano, file = "data/volcano_plot.rda", compress = "xz")
#tools::checkRdaFiles("data/volcano_plot.rda")
#tools::resaveRdaFiles("data/volcano_plot.rda")  ## why isn't this the default in 'save'?
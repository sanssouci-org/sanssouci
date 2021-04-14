library("sansSouci")

# create a temporary directory to hold the desired files
pv <- packageVersion("sansSouci")
tag <- sprintf("%s_%s", "sansSouci", pv)
dirName <- paste("python-vs-R", tag, sep = "_")
path <- file.path(tempdir(), dirName)
dir.create(path)

writeFUN <- function(x, filename, outPath = path) {
    pathname <- file.path(outPath, paste0(filename, ".csv"))
    write.table(x, file = pathname, 
                sep = ",", col.names = FALSE, row.names = FALSE)
    return(pathname)
}

data(fMRI_localizer, package = "sansSouci.data")

# only keep a small subset of the data set
m <- 543
X <- head(fMRI_localizer, m) 
fX <- writeFUN(X, "X")

# map colmun names to 0/1
categ <- ifelse(colnames(fMRI_localizer) == "left", 0, 1) 
fC <- writeFUN(categ, "categ")


# get permutation p-values
B <- 123
set.seed(0xBEEF)
p0 <- sansSouci:::get_perm_p(X, categ, B)
fperm <- writeFUN(p0, "perm_p")

# get associated pivotal statistic
pivStat <- sansSouci:::get_pivotal_stat(p0)
fpiv <- writeFUN(pivStat, "pivotal_stat")

if (FALSE) {  ## sanity check: can we read these data correctly?
    p0_bak <- as.matrix(read.table(fperm, sep = ","))
    dimnames(p0_bak) <- NULL
    stopifnot(all.equal(p0, p0_bak))
    
    pivStat_bak <- read.table(fpiv, sep = ",")[, 1]
    stopifnot(all.equal(pivStat, pivStat_bak))
}


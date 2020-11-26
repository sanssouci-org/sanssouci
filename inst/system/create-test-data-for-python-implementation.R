data(fMRI_localizer, package = "sansSouci.data")

# only keep a small subset of the data set
m <- 543
X <- head(fMRI_localizer, m) 

# map colmun names to 0/1
categ <- ifelse(colnames(fMRI_localizer) == "left", 0, 1) 

pv <- packageVersion("sansSouci")
tag <- sprintf("%s_%s", "sansSouci", pv)

# get permutation p-values
B <- 123
set.seed(0xBEEF)
p0 <- sansSouci:::get_perm_p(X, categ, B)
filename <- paste0("perm_p_", tag, ".csv")
write.table(p0, file = filename, sep = ",", col.names = FALSE, row.names = FALSE)

# get associated pivotal statistic
pivStat <- sansSouci:::get_pivotal_stat(p0)
filename <- paste0("pivotal_stat_", tag, ".csv")
write.table(pivStat, file = filename, sep = ",")

if (FALSE) {  ## sanity check: can we read these data correctly?
    p0_bak <- as.matrix(read.table("perm_p.csv", sep = ","))
    dimnames(p0_bak) <- NULL
    stopifnot(all.equal(p0, p0_bak))
    
    pivStat_bak <- read.table("pivotal_stat.csv", sep = ",")[, 1]
    stopifnot(all.equal(pivStat, pivStat_bak))
}


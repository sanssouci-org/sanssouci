library("sanssouci")

# create a temporary directory to hold the desired files
pv <- packageVersion("sansSouci")
tag <- sprintf("%s_%s", "sansSouci", pv)
dir_name <- paste("python-vs-R", tag, sep = "_")
path <- file.path(tempdir(), dir_name)
dir.create(path)

write_fun <- function(x, filename, out_path = path) {
  pathname <- file.path(out_path, paste0(filename, ".csv"))
  write.table(x,
    file = pathname,
    sep = ",", col.names = FALSE, row.names = FALSE
  )
  return(pathname)
}

data(fMRI_localizer, package = "sanssouci.data")

# only keep a small subset of the data set
m <- 543
X <- head(fMRI_localizer, m)
FX <- write_fun(X, "X")

# map colmun names to 0/1
categ <- ifelse(colnames(fMRI_localizer) == "left", 0, 1)
FC <- write_fun(categ, "categ")


# get permutation p-values
B <- 123
set.seed(0xBEEF)
p0 <- sansSouci:::get_perm_p(X, categ, B)
fperm <- write_fun(p0, "perm_p")

# get associated pivotal statistic
piv_stat <- sansSouci:::get_pivotal_stat(p0)
fpiv <- write_fun(piv_stat, "pivotal_stat")

if (FALSE) { ## sanity check: can we read these data correctly?
  p0_bak <- as.matrix(read.table(fperm, sep = ","))
  dimnames(p0_bak) <- NULL
  stopifnot(all.equal(p0, p0_bak))

  piv_stat_bak <- read.table(fpiv, sep = ",")[, 1]
  stopifnot(all.equal(piv_stat, piv_stat_bak))
}

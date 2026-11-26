data_sets <- c("chiaretti", "golub", "hedenfalk", "loi", "rosenwald")
data_sets <- data_sets[-4] ## Dropping Loi: too weak signal

get_stats <- function(data_set) {
  filename <- sprintf("%s.rds", data_set)
  pathname <- file.path("local/data", filename)
  dat <- readRDS(pathname)
  rm(pathname)

  cn <- colnames(dat)[1]
  categ <- ifelse(colnames(dat) == cn, 1, 0) # map to 0/1
  res <- row_welch_tests(dat, categ)
  mres <- Reduce(cbind, res)
  colnames(mres) <- names(res)

  m <- nrow(mres)
  gid <- rownames(mres)
  if (is.null(gid)) {
    patt <- sprintf("gene_%%0%ds", floor(log10(m)) + 1)
    gid <- sprintf(patt, 1:m)
  }
  p <- mres[, "p.value"]
  adjp <- p.adjust(p, method = "BH")
  dat <- data.frame(
    id = gid,
    dataSet = data_set,
    meanDiff = round(mres[, "meanDiff"], 3),
    p.value = p,
    adjp.BH = round(adjp, 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  dat
}

df_list <- lapply(data_sets, get_stats)
volcano <- Reduce(rbind, df_list)
save(volcano, file = "data/volcano_plot.rda", compress = "xz")

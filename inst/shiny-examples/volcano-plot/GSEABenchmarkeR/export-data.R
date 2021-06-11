library(GSEABenchmarkeR)

# change work directory to create Example data set (work also with deployed version)
orwd <- getwd()
newwd <- dir(".", "GSEABenchmarkeR", all.files = TRUE, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
setwd(paste0(getwd(),'/', newwd))

if (!dir.exists("express-data-set")){
  dir.create("express-data-set")
}
if (!dir.exists("gene-set")){
  dir.create("gene-set")
}

nameGSE <- c("GSE14762", "GSE15471", "GSE18842", "GSE19728", "GSE5281_EC",
             "GSE23878", "GSE7305", "GSE3467", "GSE9476", "GSE38666_epithelia")
data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData,"geo2kegg")
for (name in nameGSE){
  print(name)
  
  rawData <- R.cache::memoizedCall(maPreproc,data[name])[[1]]
  
  saveRDS(rawData, file=paste("express-data-set/",name,".RDS", sep = ""))
  print("ok")
}


go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
                               org = "hsa", db = "go", onto = "BP", mode = "GO.db")

cleanGo.GS <- function(go.gs){
  for(i in names(go.gs)){
    if(length(go.gs[[i]])<10){
      go.gs[i] <- NULL
    }
  }
  return(go.gs)
  
}

go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs) # our func

saveRDS(go.gs, file="gene-set/go.gs.RDS")

print("go.gs done")

setwd(orwd)
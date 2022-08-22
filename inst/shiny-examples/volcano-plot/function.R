namedGeo2kegg <- function(geo2kegg){
  namesData <- list()
  for (i in 1:length(geo2kegg)){
    x <- geo2kegg[[i]]
    ed <- Biobase::experimentData(x)
    if(!any(ed@name == c('GSE781','GSE32676'))){
      disease <- ed@other$disease
      namesData[sprintf("%s (%s)", disease, ed@name)] = ed@name
    }
  }
  return(namesData)
}

cleanGo.GS <- function(go.gs){
  for(i in names(go.gs)){
    incProgress(1/length(go.gs))
    if(length(go.gs[[i]])<6){
      go.gs[i] <- NULL
    }
  }
  return(go.gs)
  
}

readCSV_sep<- function(file, ... ){
  df <- read.csv(file, sep = c(","), ...)
  if (dim(df)[2] > 1) {
    return(df)
  }
  df <- read.csv(file, sep = c(";"), ...)
  if (dim(df)[2] > 1) {
    return(df)
  } 
  df <- read.csv(file,  sep = c("\t"), ...)
  if (dim(df)[2] > 1) {
    return(df)
  } 
  stop("Only ',', ';', '\\t' separators are allowed. Please adapte your csv file.")  
}


exampleData <- function(){
  
  #### expression matrix
  data <- list()
  data$matrix <- expr_ALL
  
  # setProgress(value = 0.4, detail = "Cleaning data ... ")
  
  categ <- colnames(data$matrix)
  data$categ <- rep(1, length(categ))
  data$categ[which(categ == "NEG")] <- 0
  
  colnames(data$matrix) <- data$categ
  
  
  # setProgress(value = 0.4, detail = "Fc and p.value matrix ... ")
  
  #### degraded matrix
  dex <- rowWelchTests(data$matrix, data$categ)
  data$degrade = data.frame("p.value"=dex[["p.value"]], "fc"=dex$estimate)
  rm(dex)
  
  #gene set matrix
  bioFun <- expr_ALL_GO
  mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
  data$biologicalFunc <- bioFun[mm, ]
  rm(bioFun)
  return(data)
}




VolcanoPlotNico <- function (X, categ, thr, p = 1, q = 1, r = 0, cex = c(0.2, 0.6), 
                             col = c("#33333333", "#FF0000", "#FF666633"), pch = 19, ylim = NULL) 
{
  if (p < 1 && q < 1) {
    warning("Filtering both on p-values and BH-adjusted p-values")
  }
  m <- nrow(X)
  dex <- rowWelchTests(X, categ)
  
  
  pval <- dex[["p.value"]]
  logp <- -log10(pval)
  fc <- dex$meanDiff
  adjp <- p.adjust(pval, method = "BH")
  y_sel <- which((adjp <= q) & (pval <= p))
  y_thr <- Inf
  if (length(y_sel) > 0) {
    y_thr <- min(logp[y_sel])
  }
  sel1 <- which(logp >= y_thr & fc >= r)
  sel2 <- which(logp >= y_thr & fc <= -r)
  sel12 <- sort(union(sel1, sel2))
  n1 <- length(sel1)
  FP1 <- maxFP(pval[sel1], thr = thr)
  TP1 <- n1 - FP1
  FDP1 <- round(FP1/max(n1, 1), 2)
  n2 <- length(sel2)
  FP2 <- maxFP(pval[sel2], thr = thr)
  TP2 <- n2 - FP2
  FDP2 <- round(FP2/max(n2, 1), 2)
  n12 <- length(sel12)
  FP12 <- maxFP(pval[sel12], thr = thr)
  TP12 <- n12 - FP12
  FDP12 <- round(FP12/max(n12, 1), 2)
  cols <- rep(col[1], m)
  cols[c(sel1, sel2)] <- col[2]
  cexs <- rep(cex[1], m)
  cexs[sel12] <- cex[2]
  xlab <- "Fold change (log scale)"
  ylab <- "p-value (-log[10] scale)"
  infty <- 100
  if (is.null(ylim)) {
    ylim <- c(0, max(logp))
  }
  lte <- "≤"
  gte <- "≥"
  title <- paste(c(sprintf("%s genes selected\nAt least %s true positives (FDP %s %s)",
                           n12, TP12, lte, FDP12)))
  txtLeft <- paste(c(sprintf("%s genes\nTP %s %s\nFDP %s %s", n2, gte,
                             TP2, lte, FDP2)))
  txtRight <- paste(c(sprintf("%s genes\nTP %s %s\nFDP %s %s", n1, gte,
                              TP1, lte, FDP1)))
  
  ggplot(data= data.frame(pvalue=logp, fc=fc, color=cols, size=cexs, stringsAsFactors=TRUE), 
         aes(x=fc, y=logp, color=as.factor(color), size = as.factor(size))) + 
    geom_point()  + 
    geom_vline(xintercept = c(-1,1)*r, colour = "grey") +
    geom_hline(yintercept = y_thr, colour = "grey") +
    scale_color_manual(values=col[1:2]) + 
    theme(legend.position='none') + 
    labs(x = xlab, y=ylab, title=title) +
    scale_size_manual(values=c(0.2,1.5)) + 
    geom_rect(aes(xmin = r, xmax = Inf, ymin = y_thr, ymax = Inf), fill = col[3], alpha = 0.002, color = NA) +
    geom_rect(aes(xmin = -r, xmax = -Inf, ymin = y_thr, ymax = Inf), fill="pink", color = NA, alpha = 0.01) + 
    annotate(x=min(fc)+0.1, y=max(logp)-1, "text", label = txtLeft) + 
    annotate(x=max(fc)-0.1, y=max(logp)-1, "text", label = txtRight)
  
}


cleanMatrix <- function(matrixFunc){
  text = ""
  boolValidation = TRUE
  color = "color:blue"
  
  for (i in 1:length(matrixFunc)){
    if(!is.numeric(matrixFunc[,i])){
      text <- "The matrix contains textual variable. We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(data = matrixFunc, text = text, boolValidation = boolValidation, color = color))
    }
  }
  
  colnam <- colnames(matrixFunc) #Control colnames(matrix)
  valueColnam <- sort(unique(colnam))
  if(length(valueColnam)==2){
    if(any(valueColnam != c("0","1"))){
      colnam[ colnam == valueColnam[1] ] <- 0
      colnam[ colnam == valueColnam[2] ] <- 1
      colnames(matrixFunc) <- colnam
      text <- paste("The colnames of your matrix does not contains 0 or 1.We consider that", valueColnam[1] ,"becomes 0 and ", valueColnam[2]," becomes 1")
      boolValidation <- TRUE
      color <- "color:orange"
    }
  }else{
    boolValidation <- FALSE
    color <- "color:red"
    text <- "The column names of your data contains more (or less) than 2 categories. Please use {0, 1} for the colnames of your matrix."
  }
  return(list(data = matrixFunc, text = text, boolValidation = boolValidation, color = color))
}

cleanBiofun <- function(biofun){
  
  text = ""
  boolValidation = TRUE
  color = "color:blue"
  
  # if biofun cotains at least one textual variable
  for (i in 1:length(biofun)){
    if(!is.numeric(biofun[,i])){
      text <- "The gene set matrix contains textual variable. We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))
    }
  }
  
  if(!all(biofun == 0 | biofun==1)){
    text <- "The gene set matrix contains other values than 0 or 1"
    boolValidation <- FALSE
    color <- "color:red"
  }
  
  return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))
}

matchMatrixBiofun <- function(geneNames, biofun){
  text = ""
  boolValidation = TRUE
  color = "color:blue"
  
  mm <- match(rownames(biofun), geneNames)
  biofun <- biofun[mm, ]
  if(all(is.na(mm))){
    text = "None of the lines of the gene set matrix correspond to the lines of the gene expression data matrix."
    boolValidation = FALSE
    color = "color:green"
  }
  
  return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))
  
}

t_linear <- function(alpha, k, m) {
  alpha * k / m;
}

calcBounds <- function(listPval, thr){
  n <- length(listPval)
  # curve_max_FP <- curveMaxFP(sort(listPval), sort(thr), "BNR2016")
  # FP <- curve_max_FP[length(curve_max_FP)]
  FP <- maxFP(sort(listPval), thr = sort(thr))
  TP <- n - FP
  FDP <- FP/max(n, 1)
  return(list(n=n, FP=FP, TP=TP, FDP=FDP))
  
}


# boundGroup <- function(df, bioFun, nameFunctions, thr){
#   table <- data.frame("Name" = c(), "# genes" = c(), "TP≥" = c(), "FDP≤"=c(), check.names = FALSE)
#   for (func in nameFunctions){
#     ids <- which(bioFun[, func] == 1)
#     listPval <- df$pval[ids]
#     bounds <- calcBounds(listPval = listPval, thr = thr)
#     table <- rbind(table, data.frame(
#       "Name" = addUrlLink(func), 
#       "# genes" = bounds$n, 
#       "TP≥" = as.integer(bounds$TP),
#       "FDP≤" = bounds$FDP,
#       check.names = FALSE))
#   }
#   table <- table[order(table["FDP≤"]),]
#   return(table)
#   
# }

# boundGroup <- function(df, bioFun, thr){
#   # create data frame
#   table <- data.frame("Name" = c(), "# genes" = c(), "TP≥" = c(), "FDP≤"=c(), check.names = FALSE)
#   # class of biofun
#   if(class(bioFun)=="list"){ 
#     nameFunctions <- names(bioFun)
#     for (func in nameFunctions){
#       incProgress(1/length(nameFunctions))
#       ids <- bioFun[[func]]
#       if (length(ids)>5){ #on ne prend que des groupes avec plus de 5 gènes
#         listPval <- df$pval[ids]
#         listPval <- listPval[!is.na(listPval)] #Calcul de la pvalue
#         if(length(listPval)>1){
#           bounds <- calcBounds(listPval = listPval, thr = thr) #calcul des bornes PH
#           table <- rbind(table, data.frame(
#             "Name" = addUrlLink(func), 
#             "# genes" = bounds$n, 
#             "TP≥" = as.integer(bounds$TP),
#             "FDP≤" = bounds$FDP,
#             check.names = FALSE))
#         }
#       } else {
#         print("length(ids)>2")
#       }
#     }
#   } else {
#     nameFunctions <- colnames(bioFun)
#     for (func in nameFunctions){
#       ids <- which(bioFun[, func] == 1)
#       listPval <- df$pval[ids]
#       bounds <- calcBounds(listPval = listPval, thr = thr)
#       table <- rbind(table, data.frame(
#         "Name" = addUrlLink(func), 
#         "# genes" = bounds$n, 
#         "TP≥" = as.integer(bounds$TP),
#         "FDP≤" = bounds$FDP,
#         check.names = FALSE))
#     }
#   }
#   table <- table[order(table["FDP≤"]),]
#   # print("boundGroup OK")
#   return(table)
#   
# }

boundGroup <- function(object){
  # create data frame
  table <- data.frame("Name" = c(), "# genes" = c(), "TP≥" = c(), "FDP≤"=c(), check.names = FALSE)
  bioFun <- object$input$biologicalFunc
  if(class(bioFun)[1]=="list"){ 
    nameFunctions <- names(bioFun)
    for (func in nameFunctions){
      incProgress(1/length(nameFunctions))
      name <- bioFun[[func]]
      ids <- which(is.element(rownames(object$input$Y),name)) 
      #le traitement des éléments de groupe non contenu dans Y est pris en compte dans is.element
      if (length(ids)>5){ #on ne prend que des groupes avec plus de 5 gènes
        bounds <- predict(object, S=ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids), 
          "TP≥" = as.integer(bounds['TP']),
          "FDP≤" = bounds['FDP'],
          check.names = FALSE))
      } 
    }
    
  } else {
    nameFunctions <- colnames(bioFun)
    for (func in nameFunctions){
      ids <- which(bioFun[, func] == 1)
      if (length(ids)>1){
        bounds <- predict(object, S=ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids), 
          "TP≥" = as.integer(bounds['TP']),
          "FDP≤" = bounds['FDP'],
          check.names = FALSE))
      }
    }
  }
  table <- table[order(table["FDP≤"]),]
  return(table)
  
}


addUrlLink <- function(name){
  if(grepl("GO:\\d+", name)){ # met ce lien spécifique si la structure GO:000000000... est dans le nom
    url <- paste("https://www.ebi.ac.uk/QuickGO/term/",str_extract_all(name, "GO:\\d+"), sep="")
    url <- paste('<a target="_blanck" href="', url, '" >',name, '</a>', sep="")
    return(url)
  } else {
    return(name)
  }
}





boundGroup2 <- function(object){
  # create data frame
  table <- data.frame("Name" = c(), "# genes" = c(), "TP≥" = c(), "FDP≤"=c(), check.names = FALSE)
  bioFun <- object$input$biologicalFunc
  if(class(bioFun)[1]=="list"){ 
    print("on passe ici")
    nameFunctions <- names(bioFun)
    for (func in nameFunctions){ 
      incProgress(1/length(nameFunctions))
      name <- bioFun[[func]]
      ids <- which(is.element(rownames(object$input$Y),name)) 
      # ids <- sapply(go.gs, function(x){which(is.element(rownames(object$input$Y),x))})
      #le traitement des éléments de groupe non contenu dans Y est pris en compte dans is.element
      #ids <- ids[sapply(ids, function(x){length(x)>5})]
      if (length(ids)>5){ #on ne prend que des groupes avec plus de 5 gènes
        bounds <- predict2(object, S=ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids), 
          "TP≥" = as.integer(bounds['TP']),
          "FDP≤" = bounds['FDP'],
          check.names = FALSE, row.names = NULL))
      } 
    }
    
  } else {
    nameFunctions <- colnames(bioFun)
    for (func in nameFunctions){
      incProgress(1/length(nameFunctions))
      ids <- which(bioFun[, func] == 1)
      if (length(ids)>1){
        bounds <- predict2(object, S=ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids), 
          "TP≥" = as.integer(bounds['TP']),
          "FDP≤" = bounds['FDP'],
          check.names = FALSE, row.names = NULL))
      }
    }
  }
  table <- table[order(table["FDP≤"]),]
  return(table)
  
}

predict2 <- function(object, S = seq_len(nHyp(object)), 
                     what = c("TP", "FDP"), all = FALSE, ...) {
  p.values <- pValues(object)
  thr <- thresholds(object)
  lab <- label(object)
  if (max(S) > nHyp(object)) {
    stop("'S' is not a subset of hypotheses")
  }
  bounds <- posthoc_bound2(p.values, S = S, thr = thr, lab = lab, what = what, all = all)
  if (!all) {
    bounds <- bounds[, "bound"]
    if (length(bounds) > 1) {
      names(bounds) <- what
    }
  }
  return(bounds)
}


posthoc_bound2 <- function (p.values, S = seq_along(p.values), thr = NULL, lab = NULL, 
                            what = c("TP", "FDP"), all = FALSE) 
{
  if (is.null(thr)) {
    stop("Argument 'thr' must be non NULL")
  }
  s <- length(S)
  idxs <- seq_len(s)
  max_FP <- rep(NA_integer_, s)
  pS <- p.values[S]
  o <- order(pS)
  sorted_p <- pS[o]
  if (length(thr) == length(p.values) && all(thr %in% c(0, 
                                                        1))) {
    max_FP <- cumsum(thr[o] == 0)
    
  }
  else {
    if (s <= 5*sqrt(length(thr))){
      max_FP <- maxFP(sorted_p, thr)
      idxs <- length(idxs)
    } else {
      max_FP <- sanssouci:::curveMaxFP(sorted_p, thr)
    }
  }
  bounds <- sanssouci:::formatBounds(max_FP, idxs = idxs, lab = lab, what = what, 
                                     all = all)
  bounds
}

































# defaultBiologicalFunc <- function(expr_ALl, expr_ALL_annotation){
#   
#   X = expr_ALL
#   categ = colnames(X)
#   categ <- rep(0, length(categ))
#   categ[colnames(X)=="NEG"] <- 1
#   cal <- calibrateJER(X, categ = categ,  B = 1e2, alpha = 0.05)
#   thr = cal$thr
#   dex <- rowWelchTests(X, categ)
#   pval <- dex[["p.value"]]
#   logp <- -log10(pval)
#   fc <- dex$meanDiff
#   adjp <- p.adjust(pval, method = "BH")
#   df = data.frame(pval, logp, fc, adjp)
#   
#   
#   A <- tibble::rownames_to_column(df, "Id")
#   
#   B = expr_ALL_annotation[c('affy_hg_u95av2','hgnc_symbol')]
#   B
#   
#   B <- dplyr::left_join(A,B, by=c("Id"="affy_hg_u95av2"))[c('Id','hgnc_symbol')]
#   rownames(B) <- B[['VALUE']]
#   # B <- B['hgnc_symbol']
#   
#   df <- tibble::rownames_to_column(df, "Id")
#   df <- dplyr::left_join(df, B, by="Id")
#   df <- df%>% 
#     rename(nameGene = hgnc_symbol)%>%
#     filter(nameGene!="")
#   rownames(df) <- df[["Id"]]
#   
#   
#   matrixFunc <- data.frame(nameGene = sort(unique(expr_ALL_annotation$hgnc_symbol)),
#                            row.names = sort(unique(expr_ALL_annotation$hgnc_symbol)))
#   matrixFunc[,c("Func1","Func2","Func3","Func4","Top4ABL1","Top19ABL1","Func7")] <- 0
#   matrixFunc[unique((df%>% filter(logp> 3) %>% filter(abs(fc)>1))[["nameGene"]]),"Func1"] <- 1
#   
#   Func2 = 712:798
#   matrixFunc[Func2, "Func2"] <- 1  
#   
#   Func3 = (df %>% filter(abs(fc)>1.4))$nameGene
#   matrixFunc[Func3, "Func3"] <- 1
#   
#   Func4 = (df %>% filter(logp > 5))$nameGene
#   matrixFunc[Func4, "Func4"] <- 1
#   
#   Func5 = c("FN1","CRK","ABL1","PPM1B","PAK4")
#   matrixFunc[Func5,"Top4ABL1"] <- 1
#   
#   Func6 = c("FN1","CRK","ABL1","PPM1B","PAK4", "TRAF6","WASL",
#             "ABI1","ABI2","CAP1","ROBO1","CBL","CRKL","INPPL1",
#             "SHC1", "HSP90AA1","DOK1","DOK2","RAD51","PXN")
#   matrixFunc[Func6,"Top19ABL1"] <- 1
#   
#   Func7 = (df %>% filter(abs(fc)<7e-4) %>% filter(logp < 2e-3))$nameGene
#   matrixFunc[Func7, "Func7"] <- 1
#   
#   return(matrixFunc)
# }

thrYaxis <- function(thr, maxlogp){
  if(maxlogp == Inf){
    maxlogp <- 16
  }
  df1 <- data.frame(num = 1:length(thr)-1, pvalue=-log10(thr))
  df2 <- data.frame(df1[c(1),])
  valeurTest <- df2[c(dim(df2)[1]),"pvalue"]
  for (i in 1:dim(df1)[1]){
    mod <- if(df1[i,"num"] < 100){ 1} else if(df1[i,"num"] < 500){ 10} else if(df1[i,"num"] < 1000){50}else{100}
    if (valeurTest - df1[i,"pvalue"] > 0.3*maxlogp/12.5 & df1[i,"num"]%%mod == 0){
      df2 <- rbind(df2, (df1[i,]))
      valeurTest <- df1[i,"pvalue"]
    }
  }
  return(df2)
}


curveMaxFP <- sanssouci:::curveMaxFP

plotMaxFP <- function(pval, thr){
  sort_pval <- sort(pval)
  curve <- curveMaxFP(sort_pval, thr, flavor="BNR2016")
  x <- 1:length(sort_pval)
  ggplot(data = data.frame(x = x, y=curve), aes(x=x, y=y)) +
    geom_line() + 
    labs(x = "Number of top features selected", y="Upper bound on the number of false positives")
}



UrlStringdbGrah <- function(vector){
  vector[2:length(vector)] <- paste0("%0d", vector[2:length(vector)])
  url <- paste("https://string-db.org/api/image/network?identifiers=", paste(vector, collapse = ""), "&species=9606", sep="")
  return(url)
}
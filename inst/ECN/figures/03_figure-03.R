## FIGURE 3

m <- nrow(Y)
q <- 0.1 # FDP budget (user-defined)
FDP <- lapply(resList, predict, what = "FDP", all = TRUE)
n_DEG <- sapply(FDP, function(x) max(which(x$bound <= q)))

library(limma)
library(edgeR)
d <- DGEList(Y)
d <- calcNormFactors(d)
Grp <- as.factor(groups)
mm <- model.matrix(~0 + Grp)

y <- voom(d, mm, plot = FALSE)

res_lm <- lmFit(y, mm)
contr <- makeContrasts(Grp1 - Grp0, levels = colnames(coef(res_lm)))
res_fit <- contrasts.fit(res_lm, contr)
res_eb <- eBayes(res_fit)
TT <- topTable(res_eb, sort.by = "none", number = Inf)

pdf("volcano-plot.pdf", width = 6, height = 6)
volcanoPlot(res, 
            fold_changes = TT$logFC, 
            p_values = TT$P.Value, 
            p = 1e-3, r = 0.5)
dev.off()

selVP <- volcanoPlot(res, 
                     fold_changes = TT$logFC, 
                     p_values = TT$P.Value, 
                     p = 1e-3, r = 0.5)
TP <- list()
for (method in names(resList)) {
    n <- n_DEG[[method]]
    TP[method] <- list(predict(resList[[method]], what = c("TP", "FDP"), S = selVP))
    #    c(Template = method, n = n, TP = TP)
}
length(selVP)
TP


## FIGURE (Appendix): limma vs Wilcoxon p-values 
df <- data.frame(wilcox = -log10(pValues(res)), limma = -log10(TT$P.Value))
p <- ggplot(df, aes(x = limma, y = wilcox)) + 
    geom_point(color = "#10101010") + 
    ggtitle("p-values (log-scale)") + 
    theme_bw() + theme(legend.position="bottom")
p
ggsave(p, file = "p-values_limma-vs-wilcoxon.pdf", width = 6, height = 4)

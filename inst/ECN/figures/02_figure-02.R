resList <- list("Adaptive Simes" = res,
                "Adaptive Simes (single step)" = res_singlestep,
                "Simes" = res_Simes,
                "ARI" = res_ARI
)

m <- nrow(Y)
q <- 0.1 # FDP budget (user-defined)
FDP <- lapply(resList, predict, what = "FDP", all = TRUE)
n_DEG <- sapply(FDP, function(x) max(which(x$bound <= q)))

TP <- NULL
for (method in names(resList)) {
    n <- n_DEG[[method]]
    TP[method] <- predict(resList[[method]], what = "TP", all = TRUE)$bound[n]
}

DEG_all <- tibble(Method = names(FDP), x = n_DEG, TP = TP, FDP = q) %>%
    pivot_longer(c(TP, FDP), names_to = "stat", values_to = "y")
conf_bounds_all <- lapply(resList, predict, all = TRUE)


## FIGURE 2
methods <- c("Adaptive Simes", "ARI")
DEG <- subset(DEG_all, Method %in% methods)
conf_bounds <- conf_bounds_all[methods]

cols <- c("black", "lightgray")

bh <- which(p.adjust(pValues(res), method = "BH") < 0.05)
xBH <- length(bh)
tpBH_res <- predict(res, what = c("TP"), all = TRUE)$bound[xBH]
fdpBH_res <- predict(res, what = c("FDP"), all = TRUE)$bound[xBH]
tpBH_resSimes <- predict(res_ARI, what = c("TP"), all = TRUE)$bound[xBH]
fdpBH_resSimes <- predict(res_ARI, what = c("FDP"), all = TRUE)$bound[xBH]

BH <- tibble(Template = "Adaptive Simes", 
             x = xBH, 
             TP = c(tpBH_res), 
             FDP = c(fdpBH_res)) %>%
    pivot_longer(c(TP, FDP), names_to = "stat", values_to = "y")

p <- plotConfCurve(conf_bounds, xmax = 2.5*max(n_DEG), cols = cols) +  
    scale_linetype_manual(values = c("solid", "solid")) +
    labs(color = "Method", linetype = "Method") +
    #    geom_vline(xintercept = n_DEG, linetype = "dotted", size = size) + 
    geom_segment(data = DEG, aes(x = x, y = -Inf, xend = x, yend = y, 
                                 color = Method, linetype = Method), 
                 size = 1) + 
    geom_segment(data = DEG, aes(x = -Inf, y = y, xend = x, yend = y, 
                                 color = Method, linetype = Method), 
                 size = 1) +
    geom_line(size = 1.2) + 
    geom_point(data = BH, aes(x = x, y = y), colour = "red", size = 2.5) +
    ggplot2::facet_wrap(~ stat, scales = "free_y") + 
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = xBH, linetype="dotted", 
               color = "red", size=0.5) 
#   + geom_hline(yintercept = q, linetype = "dashed", size = size) 
p
ggsave(p, file = "conf-curve.pdf", width = 6, height = 4)


## FIGURE 2 bis (APPENDIX)
conf_bounds <- conf_bounds_all
DEG <- DEG_all

cols <- rep(c("black", "lightgray"), each = 2)
ltys <- rep(c("solid", "dashed"), times = 2)

p <- plotConfCurve(conf_bounds, xmax = 2.5*max(n_DEG), cols = cols) +  
    scale_linetype_manual(values = ltys) +
    labs(color = "Method", linetype = "Method") +
    geom_line(size = 1.2) + 
    ggplot2::facet_wrap(~ stat, scales = "free_y") +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) 
p
ggsave(p, file = "conf-curve-annexe.pdf", width = 6, height = 4)


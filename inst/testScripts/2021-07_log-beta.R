K <- m # for m = 9038

t_beta(res_Beta$output$lambda, m-10:0, m)
t_beta_log(res_Beta_log$output$lambda, m-10:0, m)

lambda <- res_Beta$output$lambda
t_beta(lambda, m-5, m)
qbeta(lambda, m-5, m + 1 - (m-5))
1-qbeta(lambda, m-5, m + 1 - (m-5), lower.tail = FALSE);


lambda_log <- res_Beta_log$output$lambda
t_beta_log(lambda_log, m-5, m)
qbeta(lambda_log, m-5, m + 1 - (m-5), log.p = TRUE)
1-qbeta(lambda_log, m-5, m + 1 - (m-5), lower.tail = FALSE, log.p = TRUE)

# TODO fMRI data! where does K break? m =20000
library("sansSouci")
data("fMRI_localizer", package = "sansSouci.data")
Y <- fMRI_localizer
groups <- ifelse(colnames(Y) == "left", 0, 1)
obj <- SansSouci(Y = Y, groups = groups)
rm(fMRI_localizer)

m <- nrow(Y)
alpha <- 0.1
res_Simes <- fit(obj, B = 0, family = "Simes", alpha = alpha) ## B=0 => no calibration!

B <- 200
res <- fit(obj, alpha = alpha, B = B, family = "Simes")

K <- 500
res_Beta <- fit(res, alpha = alpha, B = B, family = "Beta", K = K)
res_Beta_log <- fit(res, B = B, alpha = alpha, family = "Beta-log", K = K)

max(abs(thresholds(res_Beta)-thresholds(res_Beta_log)))

resList <- list("Simes" = res_Simes,
                "Linear" = res,
                "Beta" = res_Beta,
                "Beta (log)" = res_Beta_log)
names(resList)[3] <- sprintf("Beta (K=%s)", K)
names(resList)[4] <- sprintf("Beta-log (K=%s)", K)
bounds <- sapply(resList, predict)
rownames(bounds) <- c("Lower bound on True Positives", "Upper bound on False Discovery Proportion")
knitr::kable(t(bounds), digits = 2)

bounds <- sapply(resList, predict, what = "TP")
knitr::kable(bounds, digits = 2, col.names = c("Lower bound on True Positives"))

## several K
K <- round(m/10)
res_Beta <- fit(res, alpha = alpha, B = B, family = "Beta", K = K)
res_Beta_log <- fit(res, B = B, alpha = alpha, family = "Beta-log", K = K)

df <- data.frame(template = sprintf("Beta (K=%s)", K), threshold = thresholds(res_Beta))

max(abs(thresholds(res_Beta)-thresholds(res_Beta_log)))
plot(thresholds(res_Beta)-thresholds(res_Beta_log))

K <- round(m/5)
res_Beta <- fit(res, alpha = alpha, B = B, family = "Beta", K = K)
res_Beta_log <- fit(res, B = B, alpha = alpha, family = "Beta-log", K = K)

max(abs(thresholds(res_Beta)-thresholds(res_Beta_log)))
plot(thresholds(res_Beta)-thresholds(res_Beta_log))

K <- round(m/2)
res_Beta <- fit(res, alpha = alpha, B = B, family = "Beta", K = K)
res_Beta_log <- fit(res, B = B, alpha = alpha, family = "Beta-log", K = K)

max(abs(thresholds(res_Beta)-thresholds(res_Beta_log)))
plot(thresholds(res_Beta)-thresholds(res_Beta_log))

resList <- list("Simes" = res_Simes,
                "Linear" = res,
                "Beta" = res_Beta,
                "Beta (log)" = res_Beta_log)
names(resList)[3] <- sprintf("Beta (K=%s)", K)
names(resList)[4] <- sprintf("Beta-log (K=%s)", K)
bounds <- sapply(resList, predict)
rownames(bounds) <- c("Lower bound on True Positives", "Upper bound on False Discovery Proportion")
knitr::kable(t(bounds), digits = 2)


## loop
res_list <- NULL
seq_K <- c(500, 5000, 10000, 25000)
for (ii in seq(along = seq_K)) {
    K <- seq_K[ii]
    res_Beta <- fit(res, alpha = alpha, B = B, family = "Beta", K = K)
    res_Beta_log <- fit(res, B = B, alpha = alpha, family = "Beta-log", K = K)
    df <- data.frame(family = res_Beta$parameters$family, 
                     K = K, k = 1:K,
                     threshold = thresholds(res_Beta))
    df_log <- data.frame(family = res_Beta_log$parameters$family, 
                         K = K, k = 1:K,
                         threshold = thresholds(res_Beta_log))
    res_list[[ii]] <- rbind(df, df_log)
}
length(res_list)
df <- Reduce(rbind, res_list)
df$template <- sprintf("%s (K=%s)", df$family, df$K)
head(df)

p_val <- data.frame(k = 1:m, "p_value" = sort(pValues(res_Simes)))

library("ggplot2")
ggplot(df, aes(x = k, y = threshold)) + 
    geom_line(aes(linetype = family)) +
    facet_wrap(facets = vars(K), labeller = label_both, scales = "free")
    
pl <- ggplot(df, aes(x = k, y = threshold)) + 
    geom_line(aes(color = as.factor(K), linetype = family)) + 
    geom_line(data = p_val, aes(x = k, y = p_value)) +
    xlim(c(0, 2000)) + 
    ylim(c(0, 0.01))
    
ggsave(pl, file = "beta-log_vs_p-value.png")

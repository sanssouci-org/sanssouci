for (cc in seq_configs) {
  config <- configs[cc, ]
  
  pi0 <- config[["pi0"]]
  SNR <- config[["SNR"]]
  SNR_FUN <- config[["SNR_FUN"]]
  prob <- config[["prob"]]
  
  simname <- sprintf("%s_m=%s_pi0=%s_SNR=%s%s_prob=%s_B=%s_nb-exp=%s",
                     ds_name, m, pi0, SNR_FUN, SNR, prob,
                     B, nb_exp)
  print(simname)
  cat(cc, "/", nrow(configs), ":", simname, "\n")
  filename <- sprintf("%s.rds", simname)
  pathname <- file.path(path, filename)
  
  t0 <- Sys.time()
  # res <- lapply(1:nb_exp, FUN = function(i) {
  res <- future.apply::future_lapply(1:nb_exp, future.seed = TRUE, FUN = function(i) {
    
    ## add some signal
    sig <- add_signal(X = X0, pi0 = pi0, SNR = SNR, SNR_FUN = SNR_FUN, prob = prob)
    
    ## define gene selections (for power estimation)
    tests <- rowTestFUN(sig$Y, sig$groups)
    p_values <- tests$p.value
    m <- length(p_values)
    rk <- rank(p_values)
    selections <- list(
      first_1 = which(rk %in% c(1:1)),
      first_10 = which(rk %in% c(1:10)),
      first_100 = which(rk %in% c(1:100)),
      first_1000 = which(rk %in% c(1:1000)),
      first_5000 = which(rk %in% c(1:5000)),
      first_10000 = which(rk %in% c(1:10000)),
      first_15000 = which(rk %in% c(1:15000)),
      BH_10 = which(p.adjust(p_values, method = "BH") < 0.10),
      BH_05 = which(p.adjust(p_values, method = "BH") < 0.05),
      p_05 = which(p_values < 0.05),
      p_01 = which(p_values < 0.01),
      H = 1:m)
    ## check JER control and estimate power
    res_i <- test_JER_control_DEG_intra(
      Y = sig$Y, groups = sig$groups, truth = sig$truth, 
      rowTestFUN = rowWelchTests, B = B, 
      alpha = alphas, selections = selections,
      verbose = TRUE)
    list(level = tibble(exp = i, res_i$level),
         power = tibble(exp = i, res_i$power))
  })
  level <- Reduce(rbind, lapply(res, "[[", "level"))
  power <- Reduce(rbind, lapply(res, "[[", "power"))
  res <- list(level = level, power = power)
  saveRDS(res, file = pathname)
  print(simname)
}

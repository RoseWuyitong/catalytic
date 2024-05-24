# summary for pbc data subsampling

want_index <- c(10:24)
wrong_index <- c()
for (ii in want_index) {
  load(paste0("epoch", ii, "result_Nfold10_oct4"))
  if (current_epoch_result == "this epoch contains variable with only one value") {
    print(ii)
    wrong_index <- c(wrong_index, ii)
  }
}

available_index <- want_index[-which(want_index %in% wrong_index)]


# look at all methods
mean_result_df <- data.frame()
sd_result_df <- data.frame()
pair_t_test <- matrix(0, nr = 3, nc = 4)
for (sub_i in 3:1) {
  log_PL_mean_result <- data.frame()
  MPLE_result <- data.frame()
  for (ii in available_index) {
    # load(paste0("epoch",ii,"result"))
    # load(paste0("epoch",ii,"result_Nfold10"))
    load(paste0("epoch", ii, "result_Nfold10_oct4"))
    temp <- (current_epoch_result$logPL_test_result)

    log_PL_mean_result <- rbind(log_PL_mean_result, temp[sub_i, ])
  }

  # replace MPLE result by available result instead of NA

  name_all_methods <- c(
    "MPLE", "CRE_CV", "CRE_p", "WME_CV",
    "WME_p", "post_mean_CV", "post_mean_p",
    "post_gaussian", "Ridge_CV", "Lasso_CV"
  )


  colnames(log_PL_mean_result) <- name_all_methods
  mean_result_df <- rbind(mean_result_df, apply(log_PL_mean_result, 2, mean, na.rm = TRUE))
  sd_result_df <- rbind(sd_result_df, apply(log_PL_mean_result, 2, sd, na.rm = TRUE))
  # print(apply(log_PL_mean_result,2,mean,na.rm =TRUE))
  # print(apply(log_PL_mean_result,2,sd,na.rm =TRUE))
  # print("-----------------------------------------------------------------------------------------------------------------------------------")
  # print(log_PL_mean_result)
  pair_t_test[4 - sub_i, ] <- round(c(
    t.test(log_PL_mean_result$CRE_CV, log_PL_mean_result$Ridge_CV, paired = T)$p.value,
    t.test(log_PL_mean_result$CRE_p, log_PL_mean_result$Ridge_CV, paired = T)$p.value,
    t.test(log_PL_mean_result$WME_CV, log_PL_mean_result$Ridge_CV, paired = T)$p.value,
    t.test(log_PL_mean_result$WME_p, log_PL_mean_result$Ridge_CV, paired = T)$p.value
  ), 3)
  # print(t.test(log_PL_mean_result$WME_CV,log_PL_mean_result$Ridge_CV,paired = T)$p.value)
  # print(t.test(log_PL_mean_result$CRE_CV,log_PL_mean_result$Ridge_CV,paired = T)$p.value)
}
colnames(mean_result_df) <- name_all_methods
mean_result_df
colnames(sd_result_df) <- name_all_methods
sd_result_df

merge_mean_sd <- mean_result_df
for (i in 1:3) {
  for (j in 1:10) {
    mm <- round(mean_result_df[i, j], 2)
    sdd <- round(sd_result_df[i, j] / sqrt(length(available_index)), 2)
    merge_mean_sd[i, j] <- paste0(mm, "(", sdd, ")")
  }
}
row.names(merge_mean_sd) <- c("I60", "I100", "I140")

noquote(t(merge_mean_sd))[, c(1, 2, 3)]





# pairwise t.test
noquote(t(merge_mean_sd))[c(2, 3, 4, 5, 9), ]
colnames(pair_t_test) <- c("CRE_CV", "CRE_p", "WME_CV", "WME_p")
row.names(pair_t_test) <- c("I60", "I80", "I100", "I120", "I140")
t(pair_t_test)

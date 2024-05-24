# subsampling with PM joint is forgot in previous simulation, run it now
rm(list = ls())

library(survival)
library(mvtnorm)
library(glmnet)
source("coxph_as.R")
library(emulator)
library(MASS)


get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}

generate_synthetic_X_pbc <- function(X_obs, M) {
  # column 7 and 18 are edma, need to resample together
  # categorical_index=c(1,3,4,5,6,7,18)
  syntheticMatrix <- matrix(0, nc = ncol(X_obs), nr = M)
  ###### for categorical
  for (i in c(1, 3, 4, 5, 6)) {
    u <- sample(X_obs[, i], size = M, replace = TRUE)
    ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
    u[ind] <- rbinom(sum(ind), 1, 0.5)
    syntheticMatrix[, i] <- u
  }
  #### for categorical edma
  index_resample_edma <- sample(1:nrow(X_obs), size = M, replace = T)
  ind_edma <- sample(c(TRUE, FALSE), M, replace = TRUE)
  syntheticMatrix[, c(7, 18)] <- X_obs[index_resample_edma, c(7, 18)]
  # sample (1,0) (0,1) and (0,0) random
  replace_index <- which(ind_edma)
  for (row_i in replace_index) {
    replace_value <- NULL
    ru <- runif(1)
    if (ru < 1 / 3) {
      replace_value <- c(1, 0)
    } else {
      if (ru > 2 / 3) {
        replace_value <- c(0, 1)
      } else {
        replace_value <- c(0, 0)
      }
    }
    syntheticMatrix[row_i, c(7, 18)] <- replace_value
  }
  ###### for numerical variable
  for (i in c(2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)) {
    u <- sample(X_obs[, i], size = M, replace = TRUE)
    sd_syn <- get_sd_stdnormal(X_obs[, i])
    ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
    u[ind] <- rnorm(sum(ind), mean = median(X_obs[, i]), sd = sd_syn)
    syntheticMatrix[, i] <- u
  }
  return(syntheticMatrix)
}



predict_log_PL_compute <- function(betahat, X_test, Y_test, status_test) {
  p <- ncol(X_test)
  term1 <- CoxPL(betahat, X_test, rep(0, length(Y_test)), Y_test, status_test)
  term2 <- CoxPL(rep(0, p), X_test, rep(0, length(Y_test)), Y_test, status_test)

  return(-2 * (term1 - term2))
}
source("Bayesian_cox_cat_prior_tau.R")


test_singularity_epoch <- function(subsample_index_whole, training_X, five_diff_ntrain) {
  cnt <- 0
  for (sub_i in 1:3) {
    sizetrain <- five_diff_ntrain[sub_i]
    X_ <- training_X[subsample_index_whole[1:sizetrain], ]
    num_unique <- sapply(1:ncol(X_), function(ii) {
      length(unique(X_[, ii]))
    })
    cnt <- cnt + (sum(num_unique == 1) > 0)
  }
  return(cnt)
}

one_replication_split <- function(X, Y, status, rep_index) {
  print(rep_index)
  set.seed(rep_index + 50)

  num_test <- 136
  test_index <- sample(1:nrow(X), num_test)
  training_X <- X[-test_index, ]
  test_X <- X[test_index, ]
  training_Y <- Y[-test_index]
  test_Y <- Y[test_index]
  training_status <- status[-test_index]
  test_status <- status[test_index]

  ##################### we have 140 train, give a index for subsampling
  num_train <- nrow(training_X)
  subsample_index_whole <- sample(1:num_train)
  ## we will look at 140,120,100,80,60 training size
  five_diff_ntrain_result <- rep(1000, 3)
  length_interval_gaussian_cat_cv_p <- matrix(0, nr = 3, nc = 3)
  logPL_test_result <- matrix(0, nr = 3, nc = 4)
  best_tau_cv_cat <- matrix(0, nr = 3, nc = 2)
  five_diff_ntrain <- c(140, 100, 60)
  if (test_singularity_epoch(subsample_index_whole, training_X, five_diff_ntrain) > 0) {
    return("this epoch contains variable with only one value")
  }
  for (ii in 1:3) {
    sizetrain <- five_diff_ntrain[ii]
    X_ <- training_X[subsample_index_whole[1:sizetrain], ]
    Y_ <- training_Y[subsample_index_whole[1:sizetrain]]
    status_ <- training_status[subsample_index_whole[1:sizetrain]]
    n <- nrow(X_)
    p <- ncol(X_)


    ############# generate synthetic data
    set.seed(ii * rep_index)
    M <- 1000
    X.syn <- generate_synthetic_X_pbc(X_, M)
    rate_exp <- sum(status_) / mean(Y_) / n
    T.syn <- rexp(M, rate = rate_exp)
    increase_order <- order(T.syn)
    X.syn <- X.syn[increase_order, ]
    T.syn <- T.syn[increase_order]
    status.syn <- rep(1, M)
    syntheticObj <- list()
    syntheticObj$M <- M
    syntheticObj$xstar <- X.syn
    syntheticObj$ystar <- T.syn
    syntheticObj$h_star <- rate_exp
    my_syndata <- list(T.syn = syntheticObj$ystar, X.syn = syntheticObj$xstar, status.syn = rep(1, M), h_star = rate_exp)
    ############## run some estimation

    # MPLE Ridge Lasso --------------------------------------------------------
    betahat_MLE <- as.numeric(coef(coxph(Surv(Y_, status_) ~ X_)))
    ridge_fit <- cv.glmnet(X_, Surv(Y_, status_), family = "cox", alpha = 0)
    betahat_ridgecv <- as.numeric(predict(ridge_fit, s = ridge_fit$lambda.min, type = "coeff")) #
    lasso_fit <- cv.glmnet(X_, Surv(Y_, status_), family = "cox", alpha = 1)
    betahat_lassocv <- as.numeric(predict(lasso_fit, s = lasso_fit$lambda.min, type = "coeff"))

    fit_prior_tau <- Bayesian_cox_cat_stan(X_, Y_, status_, X.syn, T.syn, status.syn)

    logPL_test_result[ii, ] <- c(
      predict_log_PL_compute(betahat_MLE, test_X, test_Y, test_status),
      predict_log_PL_compute(summary(fit_prior_tau)$summary[2:(p + 1), 1], test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_ridgecv, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_lassocv, test_X, test_Y, test_status)
    )
  }
  return(list(logPL_test_result = logPL_test_result, tau_MCMC = summary(fit_prior_tau)$summary[1, ]))
}


set.seed(111)
mydata <- pbc[complete.cases(pbc), ]
reshuffle <- sample(1:nrow(mydata), nrow(mydata))
mydata <- mydata[reshuffle, ]
X <- mydata[, 4:20]
X$trt <- ifelse(X$trt == 1, 1, 0)
X$sex <- ifelse(X$sex == "f", 1, 0)
X$edema2 <- ifelse(X$edema == 0.5, 1, 0)
X$edema <- ifelse(X$edema == 1, 1, 0)
# X$bili=log(X$bili)
# X$albumin=log(X$albumin)
# X$copper=log(X$copper)
Y <- mydata$time
Y <- Y / 365
status <- (mydata$status == 2) * 1

X <- as.matrix(X)
categorical_index <- c(1, 3, 4, 5, 6, 7, 18)
mean_seq <- apply(X[, -categorical_index], 2, mean)
sd_seq <- apply(X[, -categorical_index], 2, sd)
for (i in 1:nrow(X)) {
  X[i, -categorical_index] <- (X[i, -categorical_index] - mean_seq) / sd_seq
}




for (inded in 10:128) {
  print(paste0("current epoch is ", inded))
  current_epoch_result <- one_replication_split(X, Y, status, rep_index = inded)
  save(current_epoch_result, file = paste0("epoch", inded, "result_onlyJoint_tryhyper_oct29"))
}

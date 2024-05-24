rm(list = ls())

library(survival)
library(mvtnorm)
library(glmnet)
library(survival)
source("~/catalytic/R/COX/method_comparison/coxph_as.R")
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()



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

five_fold_cv_tau <- function(tau_seq_penalty, tau_seq_mix, X_, Y_, status_, my_syndata) {
  entry <- rep(0, nrow(X_))
  N_fold <- 10
  p <- ncol(X_)
  n <- nrow(X_)
  rep_fx_CV <- function(parallel_i) { # for CV we parallel tuning parameter at each core
    tau_penalty_i <- tau_seq_penalty[parallel_i]
    goodness_penalty_i <- 0

    tau_mix_i <- tau_seq_mix[parallel_i]
    goodness_mix_i <- 0

    five_cv_index <- lapply(1:N_fold, function(ii) {
      ((ii - 1) * n / N_fold + 1):(ii * n / N_fold)
    })
    for (i in 1:N_fold) {
      betahat_penalty_without_Bi <- betahat_pencat_function(tau_penalty_i, X_[-five_cv_index[[i]], ], Y_[-five_cv_index[[i]]], status_[-five_cv_index[[i]]], my_syndata)
      # print(betahat_penalty_without_Bi)
      term1 <- CoxPL(betahat_penalty_without_Bi, X_, entry, Y_, status_)
      term2 <- CoxPL(betahat_penalty_without_Bi, X_[-five_cv_index[[i]], ], entry[-five_cv_index[[i]]], Y_[-five_cv_index[[i]]], status_[-five_cv_index[[i]]])
      goodness_penalty_i <- goodness_penalty_i + term1 - term2



      betahat_mix_without_Bi <- betahat_mixcat_function(tau_mix_i, X_[-five_cv_index[[i]], ], Y_[-five_cv_index[[i]]], status_[-five_cv_index[[i]]], my_syndata)
      term1 <- CoxPL(betahat_mix_without_Bi, X_, entry, Y_, status_)
      term2 <- CoxPL(betahat_mix_without_Bi, X_[-five_cv_index[[i]], ], entry[-five_cv_index[[i]]], Y_[-five_cv_index[[i]]], status_[-five_cv_index[[i]]])
      goodness_mix_i <- goodness_mix_i + term1 - term2
    }
    return(c(goodness_penalty_i / N_fold, goodness_mix_i / N_fold))
  }
  NUM_TUNING <- length(tau_seq_penalty)

  all_CV_score_list <- mclapply(1:NUM_TUNING, rep_fx_CV, mc.cores = numCores)
  all_cv_score <- t(sapply(1:NUM_TUNING, function(inde) {
    all_CV_score_list[[inde]]
  }))
  best_tau_penalty_cv <- tau_seq_penalty[which.max(all_cv_score[, 1])]
  # print(all_CV_score_list)
  best_tau_mix_cv <- tau_seq_mix[which.max(all_cv_score[, 2])]
  return(list(
    best_tau_penalty = best_tau_penalty_cv, cv_score = all_cv_score,
    best_tau_mix = best_tau_mix_cv
  ))
}
betahat_pencat_function <- function(tau_this, X, Y, status, my_syndata) {
  p <- ncol(X)
  have_error22 <- FALSE
  tryCatch(
    {
      betahat_penalty <- newton_cox_syn_dec30_h_daga_const(tol = 1e-5, MAX_ITER = 25, initial_values = rep(0, p), Y, status, X, tau = tau_this, my_syndata)
    },
    error = function(e) {
      have_error22 <<- TRUE
    },
    warning = function(w) {
      have_error22 <<- TRUE
    }
  )
  if (have_error22 == TRUE) {
    have_error22 <- FALSE
    tryCatch(
      {
        betahat_penalty <- newton_cox_syn_dec30_h_daga_const(tol = 1e-5, MAX_ITER = 25, initial_values = rep(1, p), Y, status, X, tau = tau_this, my_syndata)
      },
      error = function(e) {
        have_error22 <<- TRUE
      },
      warning = function(w) {
        have_error22 <<- TRUE
      }
    )
  }
  if (have_error22 == TRUE) {
    have_error22 <- FALSE
    tryCatch(
      {
        betahat_penalty <- newton_cox_syn_dec30_h_daga_const(tol = 1e-5, MAX_ITER = 25, initial_values = rnorm(p), Y, status, X, tau = tau_this, my_syndata)
      },
      error = function(e) {
        have_error22 <<- TRUE
      },
      warning = function(w) {
        have_error22 <<- TRUE
      }
    )
  }
  if (have_error22 == TRUE) {
    return(rep(5000, ncol(X)))
  }
  return(betahat_penalty)
}

betahat_mixcat_function <- function(tau_this, X, Y, status, my_syndata) {
  M <- nrow(my_syndata$X.syn)
  X_merge <- rbind(X, my_syndata$X.syn)
  Y_merge <- c(Y, my_syndata$T.syn)
  status_merge <- c(status, my_syndata$status.syn)
  wt <- c(rep(1, nrow(X)), rep(tau_this / M, M))

  have_error22 <- FALSE
  tryCatch(
    {
      betahat_mix <- as.numeric(coef(coxph(Surv(Y_merge, status_merge) ~ X_merge, weights = wt, ties = "breslow")))
    },
    error = function(e) {
      have_error22 <<- TRUE
    },
    warning = function(w) {
      have_error22 <<- TRUE
    }
  )

  if (have_error22 == TRUE) {
    return(rep(100, ncol(X)))
  }
  return(betahat_mix)
}


predict_log_PL_compute <- function(betahat, X_test, Y_test, status_test) {
  p <- ncol(X_test)
  term1 <- CoxPL(betahat, X_test, rep(0, length(Y_test)), Y_test, status_test)
  term2 <- CoxPL(rep(0, p), X_test, rep(0, length(Y_test)), Y_test, status_test)

  return(-2 * (term1 - term2))
}

library("rstan")
setting.interval <- function(y, delta, s, J) {
  n <- length(y)

  smax <- max(s)

  case0 <- which(delta == 0)
  case1 <- which(delta == 1)

  case0yleq <- which(delta == 0 & y <= smax)
  case0ygeq <- which(delta == 0 & y > smax)
  case1yleq <- which(delta == 1 & y <= smax)
  case1ygeq <- which(delta == 1 & y > smax)


  ind.d <- ind.r <- matrix(0, n, J)

  for (i in case1yleq) {
    d.mat.ind <- min(which(s - y[i] >= 0))
    ind.d[i, d.mat.ind] <- 1
    ind.r[i, 1:d.mat.ind] <- 1
  }

  for (i in case0yleq) {
    cen.j <- min(which(s - y[i] >= 0))
    ind.r[i, 1:cen.j] <- 1
  }

  if (length(union(case1ygeq, case0ygeq)) > 0) {
    ind.r[union(case1ygeq, case0ygeq), ] <- 1
  }

  ind.r_d <- ind.r - ind.d


  d <- colSums(ind.d)

  list(
    ind.r = ind.r,
    ind.d = ind.d,
    d = d,
    ind.r_d = ind.r_d
  )
}

MODEL_cat <- stan_model("bayesian_cox_cat.stan", auto_write = rstan_options("auto_write" = TRUE))

Bayesian_cox_cat_stan <- function(MODEL_cat, X, Y, status, X.syn, T.syn, status.syn, tau = NULL, h_daga = NULL) {
  # Extract dimensions
  p <- ncol(X)
  n <- nrow(X)
  M <- nrow(X.syn)

  # Set default values
  if (is.null(tau)) {
    tau <- p
  }

  if (is.null(h_daga)) {
    h_daga <- function(t) {
      sum(status) / mean(Y) / n
    }
  }

  # Hyperparameters for gamma process prior construction
  eta0 <- 1
  kappa0 <- 1

  weibull_fit <- survreg(Surv(Y, status) ~ 1, dist = "weibull")
  shape_fit <- 1 / weibull_fit$scale
  scale_fit <- exp(coef(weibull_fit))
  c0 <- 2
  s <- c(sort(Y[status == 1]), 2 * max(Y) - max(Y[-which(Y == max(Y))]))
  s <- unique(s)
  J <- length(s)

  # Calculate intervals
  intv <- setting.interval(Y, status, s, J)

  # Compute H.star and alpha0
  H.star <- alpha0 <- numeric(J)
  for (j in 1:J) {
    H.star[j] <- (s[j] / scale_fit)^shape_fit # eta0 * s[j] ^ kappa0
    alpha0[j] <- c0 * H.star[j]
  }

  hPriorSh <- diff(c(0, alpha0))

  # Compute H_daga_Y_star
  H_daga_Y_star <- sapply(1:M, function(index) {
    integrate(Vectorize(h_daga), lower = 0, upper = T.syn[index])$value
  })

  # For grouped likelihood
  R_tilde_minus_D_tilde <- intv$ind.r_d
  D_tilde <- intv$ind.d

  # Stan model and sampling

  fit <- sampling(
    MODEL_cat,
    data = list(
      J = J,
      hPriorSh = hPriorSh,
      c0 = c0,
      eta0 = eta0,
      kappa0 = kappa0,
      n = n,
      p = p,
      X = X,
      M = M,
      X_star = X.syn,
      Y_star = T.syn,
      status_star = status.syn,
      H_daga_Y_star = H_daga_Y_star,
      tau_downweight = tau,
      R_tilde_minus_D_tilde = R_tilde_minus_D_tilde,
      D_tilde = D_tilde
    ),
    chains = 1
  )

  return(fit)
}



MODEL_gaussian <- stan_model("bayesian_cox_gaussian.stan", auto_write = rstan_options("auto_write" = TRUE))

Bayesian_cox_gaussian_stan <- function(MODEL_gaussian, X, Y, status, lam = NULL, h_daga = NULL) {
  # Extract dimensions
  p <- ncol(X)
  n <- nrow(X)


  # Set default values
  if (is.null(lam)) {
    lam <- p
  }

  if (is.null(h_daga)) {
    h_daga <- function(t) {
      sum(status) / mean(Y) / n
    }
  }

  # Hyperparameters for gamma process prior construction
  eta0 <- 1
  kappa0 <- 1

  weibull_fit <- survreg(Surv(Y, status) ~ 1, dist = "weibull")
  shape_fit <- 1 / weibull_fit$scale
  scale_fit <- exp(coef(weibull_fit))
  c0 <- 2
  s <- c(sort(Y[status == 1]), 2 * max(Y) - max(Y[-which(Y == max(Y))]))
  s <- unique(s)
  J <- length(s)

  # Calculate intervals
  intv <- setting.interval(Y, status, s, J)

  # Compute H.star and alpha0
  H.star <- alpha0 <- numeric(J)
  for (j in 1:J) {
    H.star[j] <- (s[j] / scale_fit)^shape_fit # eta0 * s[j] ^ kappa0
    alpha0[j] <- c0 * H.star[j]
  }

  hPriorSh <- diff(c(0, alpha0))



  # For grouped likelihood
  R_tilde_minus_D_tilde <- intv$ind.r_d
  D_tilde <- intv$ind.d

  # Stan model and sampling

  fit <- sampling(
    MODEL_gaussian,
    data = list(
      J = J,
      hPriorSh = hPriorSh,
      c0 = c0,
      eta0 = eta0,
      kappa0 = kappa0,
      n = n,
      p = p,
      X = X,
      lambda = lam,
      R_tilde_minus_D_tilde = R_tilde_minus_D_tilde,
      D_tilde = D_tilde
    ),
    chains = 1
  )

  return(fit)
}


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
  ## we will look at 140,100,60 training size
  five_diff_ntrain_result <- rep(1000, 3)
  length_interval_gaussian_cat_cv_p <- matrix(0, nr = 3, nc = 3)
  logPL_test_result <- matrix(0, nr = 3, nc = 10)
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

    tau_seq_penalty <- 10^seq(-1, log10(nrow(X_)), length.out = NUM_TUNING)
    tau_seq_mix <- 10^seq(-1, log10(nrow(X_)), length.out = NUM_TUNING)

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
    # ridge_fit$lambda.min is best lambda for ridge, will plug in for gaussian prior
    # need to use lambda.min from standardize=false


    # point estimation with synthetic data ------------------------------------
    cv_result_list <- five_fold_cv_tau(tau_seq_penalty, tau_seq_mix, X_, Y_, status_, my_syndata)
    best_tau_penalty_cv <- cv_result_list$best_tau_penalty
    best_tau_mix_cv <- cv_result_list$best_tau_mix
    betahat_tau_cv_mix <- betahat_mixcat_function(best_tau_mix_cv, X_, Y_, status_, my_syndata) #
    betahat_tau_cv <- betahat_pencat_function(best_tau_penalty_cv, X_, Y_, status_, my_syndata) #
    betahat_tau_p <- betahat_pencat_function(p, X_, Y_, status_, my_syndata)
    betahat_tau_p_mix <- betahat_mixcat_function(p, X_, Y_, status_, my_syndata) #

    best_tau_cv_cat[ii, ] <- c(best_tau_penalty_cv, best_tau_mix_cv)
    # Baye: gaussian prior, cat_cv prior, cat_p prior
    infoX <- scale(X_)
    meann <- attr(infoX, "scaled:center")
    sdd <- attr(infoX, "scaled:scale")
    fit_gaussian_prior <- Bayesian_cox_gaussian_stan(MODEL_gaussian, scale(X_), Y_, status_, lam = n * ridge_fit$lambda.min)
    fit_cat_cv <- Bayesian_cox_cat_stan(MODEL_cat, X_, Y_, status_, X.syn, T.syn, status.syn, tau = best_tau_penalty_cv)
    fit_cat_p <- Bayesian_cox_cat_stan(MODEL_cat, X_, Y_, status_, X.syn, T.syn, status.syn, tau = p)

    length_interval_gaussian_cat_cv_p[ii, ] <- c(
      mean(summary(fit_gaussian_prior)$summary[1:p, "97.5%"] / sdd - summary(fit_gaussian_prior)$summary[1:p, "2.5%"] / sdd),
      mean(summary(fit_cat_cv)$summary[1:p, "97.5%"] - summary(fit_cat_cv)$summary[1:p, "2.5%"]),
      mean(summary(fit_cat_p)$summary[1:p, "97.5%"] - summary(fit_cat_p)$summary[1:p, "2.5%"])
    )

    logPL_test_result[ii, ] <- c(
      predict_log_PL_compute(betahat_MLE, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_tau_cv, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_tau_p, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_tau_cv_mix, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_tau_p_mix, test_X, test_Y, test_status),
      predict_log_PL_compute(summary(fit_cat_cv)$summary[1:p, 1], test_X, test_Y, test_status),
      predict_log_PL_compute(summary(fit_cat_p)$summary[1:p, 1], test_X, test_Y, test_status),
      predict_log_PL_compute(summary(fit_gaussian_prior)$summary[1:p, 1] / sdd, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_ridgecv, test_X, test_Y, test_status),
      predict_log_PL_compute(betahat_lassocv, test_X, test_Y, test_status)
    )
  }
  return(list(
    logPL_test_result = logPL_test_result, length_interval_gaussian_cat_cv_p = length_interval_gaussian_cat_cv_p,
    best_tau_cv_cat = best_tau_cv_cat
  ))
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




NUM_TUNING <- 100

N_epoch <- 70


for (inded in 10:128) {
  print(paste0("current epoch is ", inded))
  current_epoch_result <- one_replication_split(X, Y, status, rep_index = inded)
  save(current_epoch_result, file = paste0("epoch", inded, "result_Nfold10_oct29"))
}

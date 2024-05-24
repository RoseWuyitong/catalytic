rm(list = ls())
library(mvtnorm)
library(glmnet)
library(survival)
source("coxph_as.R")
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()
get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}
generate_synthetic_X <- function(X_obs, M) {
  apply(X_obs, 2, function(v) {
    u <- sample(v, size = M, replace = TRUE)
    if (length(unique(v)) <= 2) {
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rbinom(sum(ind), 1, 0.5)
    } else {
      sd_syn <- get_sd_stdnormal(v)
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rnorm(sum(ind), mean = median(v), sd = sd_syn)
    }
    u
  })
}

N_rep <- 100
M <- 1000
n_seq <- c(100)

NUM_TUNING <- 100

h0_indication <- "const"
rho <- 0.5
censor_rate <- 0.1
p <- 40

five_fold_cv_tau <- function(tau_seq_penalty, tau_seq_mix, X, Y, status, my_syndata, censor_rate, rho) {
  entry <- rep(0, nrow(X))
  N_fold <- 10
  p <- ncol(X)
  n <- nrow(X)
  rep_fx_CV <- function(parallel_i) { # for CV we parallel tuning parameter at each core
    tau_penalty_i <- tau_seq_penalty[parallel_i]
    goodness_penalty_i <- 0

    tau_mix_i <- tau_seq_mix[parallel_i]
    goodness_mix_i <- 0

    five_cv_index <- lapply(1:N_fold, function(i) {
      ((i - 1) * n / N_fold + 1):(i * n / N_fold)
    })
    for (i in 1:N_fold) {
      have_error <- FALSE
      tryCatch(
        {
          betahat_penalty_without_Bi <- betahat_pencat_function(tau_penalty_i, X[-five_cv_index[[i]], ], Y[-five_cv_index[[i]]], status[-five_cv_index[[i]]], my_syndata)
        },
        error = function(e) {
          have_error <<- TRUE
        },
        warning = function(w) {
          have_error <<- TRUE
        }
      )

      if (have_error == TRUE) {
        goodness_penalty_i <- 0 - 1000000000
      } else {
        term1 <- CoxPL(betahat_penalty_without_Bi, X, entry, Y, status)
        term2 <- CoxPL(betahat_penalty_without_Bi, X[-five_cv_index[[i]], ], entry[-five_cv_index[[i]]], Y[-five_cv_index[[i]]], status[-five_cv_index[[i]]])
        goodness_penalty_i <- goodness_penalty_i + term1 - term2
      }


      have_error <- FALSE
      tryCatch(
        {
          betahat_mix_without_Bi <- betahat_mixcat_function(tau_mix_i, X[-five_cv_index[[i]], ], Y[-five_cv_index[[i]]], status[-five_cv_index[[i]]], my_syndata)
        },
        error = function(e) {
          have_error <<- TRUE
        },
        warning = function(w) {
          have_error <<- TRUE
        }
      )

      if (have_error == TRUE) {
        goodness_mix_i <- 0 - 1000000000
      } else {
        term1 <- CoxPL(betahat_mix_without_Bi, X, entry, Y, status)
        term2 <- CoxPL(betahat_mix_without_Bi, X[-five_cv_index[[i]], ], entry[-five_cv_index[[i]]], Y[-five_cv_index[[i]]], status[-five_cv_index[[i]]])
        goodness_mix_i <- goodness_mix_i + term1 - term2
      }
    }
    return(c(goodness_penalty_i / N_fold, goodness_mix_i / N_fold))
  }
  NUM_TUNING <- length(tau_seq_penalty)

  all_CV_score_list <- mclapply(1:NUM_TUNING, rep_fx_CV, mc.cores = numCores)
  all_cv_score <- t(sapply(1:NUM_TUNING, function(inde) {
    all_CV_score_list[[inde]]
  }))
  best_tau_penalty_cv <- tau_seq_penalty[which.max(all_cv_score[, 1])]
  best_tau_mix_cv <- tau_seq_mix[which.max(all_cv_score[, 2])]
  return(list(
    best_tau_penalty = best_tau_penalty_cv, cv_score = all_cv_score,
    best_tau_mix = best_tau_mix_cv
  ))
}
betahat_pencat_function <- function(tau_this, X, Y, status, my_syndata) {
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

betahat_MPLE_store <- matrix(0, nr = N_rep, nc = p)
betahat_ridge_cv_store <- matrix(0, nr = N_rep, nc = p)
betahat_lasso_cv_store <- matrix(0, nr = N_rep, nc = p)
betahat_pencat_cv_store <- matrix(0, nr = N_rep, nc = p)
betahat_mixcat_cv_store <- matrix(0, nr = N_rep, nc = p)
betahat_pencat_taup_store <- matrix(0, nr = N_rep, nc = p)




for (n in n_seq) {
  tau_seq_penalty <- 10^seq(-2, log10(n), length.out = NUM_TUNING)
  tau_seq_mix <- 10^seq(-4, log10(n), length.out = NUM_TUNING)
  FILENAME <- paste0("p", p, "censor", censor_rate * 10, "n", n)
  print(FILENAME)
  all_result_squareloss_matrix <- matrix(0, nrow = N_rep, nc = 11)
  for (rep_index in 1:N_rep) {
    set.seed(50 + rep_index)

    X_obs <- matrix(rnorm(n * p), n, p)
    # X_obs<-scale(X_obs)
    X_obs[, 1] <- rbinom(n, 1, 0.9)

    X_obs[, 3] <- rchisq(n, 1)
    X_obs[, 2] <- rchisq(n, 4)

    X <- X_obs
    beta <- c(4, -4, 3, -3, 1, -1, 1, -1, rep(1, p - 8)) / sqrt(p)

    Xbeta <- as.numeric(X %*% beta)
    TT <- sapply(1:n, function(i) {
      rexp(1, 0.5 * exp(Xbeta[i]))
    })
    if (censor_rate == 0.1) {
      pre_value <- quantile(TT, 0.95)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.2) {
      pre_value <- quantile(TT, 0.90)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.4) {
      pre_value <- quantile(TT, 0.70)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    Y_status <- t(sapply(1:n, function(aa) {
      c(min(C[aa], TT[aa]), (C[aa] > TT[aa]))
    }))
    Y <- Y_status[, 1]
    status <- Y_status[, 2]


    have_warning_MPLE <- FALSE
    betahat_MLE <- tryCatch(
      {
        as.numeric(coef(coxph(Surv(Y, status) ~ X)))
      },
      error = function(e) {
        have_warning_MPLE <<- TRUE
      },
      warning = function(w) {
        have_warning_MPLE <<- TRUE
        rep(5000, p)
      }
    )
    ###### ridge
    ridge_fit <- cv.glmnet(X, Surv(Y, status), family = "cox", alpha = 0, nfolds = 10)
    betahat_ridgecv <- as.numeric(predict(ridge_fit, s = ridge_fit$lambda.min, type = "coeff"))
    ###### lasso
    lasso_fit <- cv.glmnet(X, Surv(Y, status), family = "cox", alpha = 1, nfolds = 10)
    betahat_lassocv <- as.numeric(predict(lasso_fit, s = lasso_fit$lambda.min, type = "coeff"))




    X.syn <- X.syn <- generate_synthetic_X(X, M)
    rate_exp <- sum(status) / mean(Y) / n
    T.syn <- rexp(M, rate = rate_exp)
    increase_order <- order(T.syn)
    X.syn <- X.syn[increase_order, ]
    T.syn <- T.syn[increase_order]
    syntheticObj <- list()
    syntheticObj$M <- M
    syntheticObj$xstar <- X.syn
    syntheticObj$ystar <- T.syn
    syntheticObj$h_star <- rate_exp
    my_syndata <- list(T.syn = syntheticObj$ystar, X.syn = syntheticObj$xstar, status.syn = rep(1, M), h_star = rate_exp)
    ## CV tau
    cv_result_list <- five_fold_cv_tau(tau_seq_penalty, tau_seq_mix, X, Y, status, my_syndata, censor_rate, rho)
    best_tau_penalty_cv <- cv_result_list$best_tau_penalty
    best_tau_mix_cv <- cv_result_list$best_tau_mix
    cv_score <- cv_result_list$cv_score[, 1]

    betahat_tau_cv_mix <- betahat_mixcat_function(best_tau_mix_cv, X, Y, status, my_syndata)
    betahat_tau_cv <- betahat_pencat_function(best_tau_penalty_cv, X, Y, status, my_syndata)
    betahat_tau_p <- betahat_pencat_function(p, X, Y, status, my_syndata)


    # what we want to return and save: square loss, all beta  , wheher MPLE break up or not
    all_result_squareloss_matrix[rep_index, ] <- c(
      norm(beta - betahat_MLE, "2")^2,
      norm(beta - betahat_ridgecv, "2")^2,
      norm(beta - betahat_lassocv, "2")^2,
      norm(beta - betahat_tau_cv, "2")^2,
      norm(beta - betahat_tau_p, "2")^2,
      norm(beta - betahat_tau_cv_mix, "2")^2,
      have_warning_MPLE, best_tau_penalty_cv, best_tau_mix_cv, p, n
    )
    betahat_pencat_taup_store[rep_index, ] <- betahat_tau_p
    betahat_pencat_cv_store[rep_index, ] <- betahat_tau_cv
    betahat_mixcat_cv_store[rep_index, ] <- betahat_tau_cv_mix
    betahat_ridge_cv_store[rep_index, ] <- betahat_ridgecv
    betahat_lasso_cv_store[rep_index, ] <- betahat_lassocv
    betahat_MPLE_store[rep_index, ] <- betahat_MLE
  }
  save(all_result_squareloss_matrix, file = paste0(FILENAME, "square_loss"))
  save(betahat_MPLE_store, file = paste0(FILENAME, "betahat_MPLE"))
  save(betahat_ridge_cv_store, file = paste0(FILENAME, "betahat_ridgecv"))
  save(betahat_lasso_cv_store, file = paste0(FILENAME, "betahat_lassocv"))
  save(betahat_pencat_cv_store, file = paste0(FILENAME, "betahat_taucv"))
  save(betahat_mixcat_cv_store, file = paste0(FILENAME, "betahat_taumixcv"))
  save(betahat_pencat_taup_store, file = paste0(FILENAME, "betahat_taup"))
}

rm(list = ls())
library(mvtnorm)
library(glmnet)
library(survival)
source("coxph_as.R")
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()
get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}
generate_synthetic_X <- function(X_obs, M) {
  apply(X_obs, 2, function(v) {
    u <- sample(v, size = M, replace = TRUE)
    if (length(unique(v)) <= 2) {
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rbinom(sum(ind), 1, 0.5)
    } else {
      sd_syn <- get_sd_stdnormal(v)
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rnorm(sum(ind), mean = median(v), sd = sd_syn)
    }
    u
  })
}

N_rep <- 100
M <- 1000
n_seq <- c(100)

NUM_TUNING <- 100

h0_indication <- "const"
rho <- 0.5
censor_rate <- 0.1
p <- 40




# start bayesian sampling -------------------------------------------------

betahat_coverage_cv_store <- rep(0, p)
betahat_coverage_taup_store <- rep(0, p)
betahat_length_interval_cv_store <- rep(0, p)
betahat_length_interval_taup_store <- rep(0, p)
betahat_baye_cv_store <- matrix(0, nr = N_rep, nc = p)
betahat_baye_taup_store <- matrix(0, nr = N_rep, nc = p)

betahat_coverage_gaussian_store <- rep(0, p)
betahat_length_interval_gaussian_store <- rep(0, p)
betahat_baye_gaussian_store <- matrix(0, nr = N_rep, nc = p)




source("Bayesian_cox_cat.R")

for (n in n_seq) {
  FILENAME <- paste0("p", p, "censor", censor_rate * 10, "n", n)
  print(FILENAME)
  load(paste0(FILENAME, "square_loss"))
  rep_f <- function(rep_index) {
    set.seed(50 + rep_index)
    print(paste0(n, "rep:", rep_index))
    X_obs <- matrix(rnorm(n * p), n, p)
    # X_obs<-scale(X_obs)
    X_obs[, 1] <- rbinom(n, 1, 0.9)

    X_obs[, 3] <- rchisq(n, 1)
    X_obs[, 2] <- rchisq(n, 4)

    X <- X_obs
    beta <- c(4, -4, 3, -3, 1, -1, 1, -1, rep(1, p - 8)) / sqrt(p)

    Xbeta <- as.numeric(X %*% beta)
    TT <- sapply(1:n, function(i) {
      rexp(1, 0.5 * exp(Xbeta[i]))
    })
    if (censor_rate == 0.1) {
      pre_value <- quantile(TT, 0.95)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.2) {
      pre_value <- quantile(TT, 0.90)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.4) {
      pre_value <- quantile(TT, 0.70)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    Y_status <- t(sapply(1:n, function(aa) {
      c(min(C[aa], TT[aa]), (C[aa] > TT[aa]))
    }))
    Y <- Y_status[, 1]
    status <- Y_status[, 2]



    X.syn <- X.syn <- generate_synthetic_X(X, M)
    rate_exp <- sum(status) / mean(Y) / n
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
    if (all_result_squareloss_matrix[rep_index, 7] == 1) {
      return(rbind(rep(2000, p), rep(2000, p), rep(2000, p), rep(1000, p), rep(1000, p), rep(1000, p), rep(1000, p), rep(1000, p), rep(1000, p)))
    } else {
      best_tau_penalty_cv <- all_result_squareloss_matrix[rep_index, 8]
      fitcv <- Bayesian_cox_cat_stan(MODEL, X, Y, status, X.syn, T.syn, status.syn, tau = best_tau_penalty_cv)
      LOWER95cv <- summary(fitcv)$summary[1:p, "2.5%"]
      UPPER95cv <- summary(fitcv)$summary[1:p, "97.5%"]

      fitp <- Bayesian_cox_cat_stan(MODEL, X, Y, status, X.syn, T.syn, status.syn, tau = p)
      LOWER95p <- summary(fitp)$summary[1:p, "2.5%"]
      UPPER95p <- summary(fitp)$summary[1:p, "97.5%"]

      infoX <- scale(X)
      meann <- attr(infoX, "scaled:center")
      sdd <- attr(infoX, "scaled:scale")
      ridge_fitstdF <- cv.glmnet(scale(X), Surv(Y, status), family = "cox", alpha = 0, standardize = F)
      fit_gaussian_prior <- fitp # Bayesian_cox_gaussian_stan(MODEL_gaussian,scale(X), Y, status,lam=n*(ridge_fitstdF$lambda.min))
      LOWER95gaussian <- summary(fit_gaussian_prior)$summary[1:p, "2.5%"] / sdd
      UPPER95gaussian <- summary(fit_gaussian_prior)$summary[1:p, "97.5%"] / sdd


      return(rbind(
        summary(fitcv)$summary[1:p, 1], summary(fitp)$summary[1:p, 1], summary(fit_gaussian_prior)$summary[1:p, 1] / sdd,
        rep(mean(beta < UPPER95cv & beta > LOWER95cv), p),
        rep(mean(beta < UPPER95p & beta > LOWER95p), p),
        rep(mean(beta < UPPER95gaussian & beta > LOWER95gaussian), p),
        UPPER95cv - LOWER95cv,
        UPPER95p - LOWER95p,
        UPPER95gaussian - LOWER95gaussian
      ))
    }
  }
  betahat_baye_list <- mclapply(1:N_rep, rep_f, mc.cores = numCores)
  for (ii in 1:N_rep) {
    betahat_baye_cv_store[ii, ] <- betahat_baye_list[[ii]][1, ]
    betahat_baye_taup_store[ii, ] <- betahat_baye_list[[ii]][2, ]
    betahat_baye_gaussian_store[ii, ] <- betahat_baye_list[[ii]][3, ]
    betahat_coverage_cv_store[ii] <- betahat_baye_list[[ii]][4, 1]
    betahat_coverage_taup_store[ii] <- betahat_baye_list[[ii]][5, 1]
    betahat_coverage_gaussian_store[ii] <- betahat_baye_list[[ii]][6, 1]
    betahat_length_interval_cv_store[ii] <- mean(betahat_baye_list[[ii]][7, ])
    betahat_length_interval_taup_store[ii] <- mean(betahat_baye_list[[ii]][8, ])
    betahat_length_interval_gaussian_store[ii] <- mean(betahat_baye_list[[ii]][9, ])
  }
  save(betahat_baye_cv_store, file = paste0(FILENAME, "bayefitmean_cv"))
  save(betahat_baye_taup_store, file = paste0(FILENAME, "bayefitmean_p"))
  save(betahat_coverage_cv_store, file = paste0(FILENAME, "bayefitcoverage_cv"))
  save(betahat_coverage_taup_store, file = paste0(FILENAME, "bayefitcoverage_p"))
  save(betahat_length_interval_cv_store, file = paste0(FILENAME, "bayefitclength_cv"))
  save(betahat_length_interval_taup_store, file = paste0(FILENAME, "bayefitlength_p"))

  save(betahat_baye_gaussian_store, file = paste0(FILENAME, "bayefitmean_gaussian"))

  save(betahat_coverage_gaussian_store, file = paste0(FILENAME, "bayefitcoverage_gaussian"))

  save(betahat_length_interval_gaussian_store, file = paste0(FILENAME, "bayefitclength_gaussian"))
}




rm(list = ls())
library(mvtnorm)
library(glmnet)
library(survival)
source("coxph_as.R")
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()
get_sd_stdnormal <- function(ddd) {
  sigma_seq <- seq(0.1, 10, by = 0.02)
  target_index <- which.min(abs(pnorm(IQR(ddd) / 2, mean = 0, sd = sigma_seq) - 0.75))
  return(sigma_seq[target_index])
}
generate_synthetic_X <- function(X_obs, M) {
  apply(X_obs, 2, function(v) {
    u <- sample(v, size = M, replace = TRUE)
    if (length(unique(v)) <= 2) {
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rbinom(sum(ind), 1, 0.5)
    } else {
      sd_syn <- get_sd_stdnormal(v)
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rnorm(sum(ind), mean = median(v), sd = sd_syn)
    }
    u
  })
}

N_rep <- 100
M <- 1000
n_seq <- c(100)

NUM_TUNING <- 100

h0_indication <- "const"
rho <- 0.5
censor_rate <- 0.1
p <- 40

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


betahat_baye_tauprior_store <- matrix(0, nr = N_rep, nc = p)
beta_coverage_tauprior_store <- rep(0, N_rep)
betahat_length_interval_tauprior_store <- rep(0, p)
betahat_coverage_gaussian_store <- rep(0, p)
betahat_length_interval_gaussian_store <- rep(0, p)
betahat_baye_gaussian_store <- matrix(0, nr = N_rep, nc = p)


source("Bayesian_cox_cat_prior_tau.R")


for (n in n_seq) {
  FILENAME <- paste0("p", p, "censor", censor_rate * 10, "n", n)
  print(FILENAME)
  load(paste0(FILENAME, "square_loss"))


  rep_f <- function(rep_index) {
    set.seed(50 + rep_index)
    print(paste0(n, "rep:", rep_index))
    X_obs <- matrix(rnorm(n * p), n, p)
    # X_obs<-scale(X_obs)
    X_obs[, 1] <- rbinom(n, 1, 0.9)

    X_obs[, 3] <- rchisq(n, 1)
    X_obs[, 2] <- rchisq(n, 4)

    X <- X_obs
    beta <- c(4, -4, 3, -3, 1, -1, 1, -1, rep(1, p - 8)) / sqrt(p)

    Xbeta <- as.numeric(X %*% beta)
    TT <- sapply(1:n, function(i) {
      rexp(1, 0.5 * exp(Xbeta[i]))
    })
    if (censor_rate == 0.1) {
      pre_value <- quantile(TT, 0.95)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.2) {
      pre_value <- quantile(TT, 0.90)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    if (censor_rate == 0.4) {
      pre_value <- quantile(TT, 0.70)

      C <- runif(n, 0, pre_value)
      censor_rates <- sum(C < TT) / n
      target_censor_rate <- censor_rate
      while (abs(target_censor_rate - censor_rates) > 0.02) {
        C <- runif(n, 0, pre_value)
        censor_rates <- sum(C < TT) / n
        # print(censor_rates)
      }
    }
    Y_status <- t(sapply(1:n, function(aa) {
      c(min(C[aa], TT[aa]), (C[aa] > TT[aa]))
    }))
    Y <- Y_status[, 1]
    status <- Y_status[, 2]



    X.syn <- X.syn <- generate_synthetic_X(X, M)
    rate_exp <- sum(status) / mean(Y) / n
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

    betahat_tau_p_mix <- betahat_mixcat_function(p, X, Y, status, my_syndata)
    # betahat_mixcat_p_store[rep_index,]=betahat_tau_p_mix
    if (all_result_squareloss_matrix[rep_index, 7] == 1) {
      return(rbind(rep(2000, p), rep(2000, p), rep(2000, p), rep(2000, p), rep(2000, p), rep(2000, p)))
    } else {
      have_error22 <- FALSE
      tryCatch(
        {
          fit_prior_tau <- Bayesian_cox_cat_stan(X, Y, status, X.syn, T.syn, status.syn)
        },
        error = function(e) {
          have_error22 <<- TRUE
        },
        warning = function(w) {
          have_error22 <<- TRUE
        }
      )
      if (have_error22 == TRUE) {
        set.seed(666)
        fit_prior_tau <- Bayesian_cox_cat_stan(X, Y, status, X.syn, T.syn, status.syn)
      }

      LOWER95 <- summary(fit_prior_tau)$summary[2:(p + 1), "2.5%"]
      UPPER95 <- summary(fit_prior_tau)$summary[2:(p + 1), "97.5%"]

      infoX <- scale(X)
      meann <- attr(infoX, "scaled:center")
      sdd <- attr(infoX, "scaled:scale")
      ridge_fitstdF <- cv.glmnet(X, Surv(Y, status), family = "cox", alpha = 0)
      fit_gaussian_prior <- Bayesian_cox_gaussian_stan(MODEL_gaussian, scale(X), Y, status, lam = n * (ridge_fitstdF$lambda.min))
      LOWER95gaussian <- summary(fit_gaussian_prior)$summary[1:p, "2.5%"] / sdd
      UPPER95gaussian <- summary(fit_gaussian_prior)$summary[1:p, "97.5%"] / sdd

      return(rbind(
        summary(fit_prior_tau)$summary[2:(p + 1), 1],
        rep(mean(beta < UPPER95 & beta > LOWER95), p), UPPER95 - LOWER95,
        summary(fit_gaussian_prior)$summary[1:p, 1] / sdd,
        rep(mean(beta < UPPER95gaussian & beta > LOWER95gaussian), p),
        UPPER95gaussian - LOWER95gaussian
      ))
    }
  }

  three_result_list <- mclapply(1:N_rep, rep_f, mc.cores = numCores)
  for (rep_index in 1:N_rep) {
    betahat_baye_tauprior_store[rep_index, ] <- three_result_list[[rep_index]][1, ]
    beta_coverage_tauprior_store[rep_index] <- three_result_list[[rep_index]][2, 1]
    betahat_length_interval_tauprior_store[rep_index] <- mean(three_result_list[[rep_index]][3, ])

    betahat_baye_gaussian_store[rep_index, ] <- three_result_list[[rep_index]][4, ]
    betahat_coverage_gaussian_store[rep_index] <- three_result_list[[rep_index]][5, 1]
    betahat_length_interval_gaussian_store[rep_index] <- mean(three_result_list[[rep_index]][6, ])
  }


  save(betahat_baye_tauprior_store, file = paste0(FILENAME, "bayefitmean_tauprior"))
  save(beta_coverage_tauprior_store, file = paste0(FILENAME, "bayefitcoverage_tauprior"))
  save(betahat_length_interval_tauprior_store, file = paste0(FILENAME, "bayefitlength_tauprior"))
  save(betahat_baye_gaussian_store, file = paste0(FILENAME, "bayefitmean_gaussian"))
  save(betahat_coverage_gaussian_store, file = paste0(FILENAME, "bayefitcoverage_gaussian"))
  save(betahat_length_interval_gaussian_store, file = paste0(FILENAME, "bayefitclength_gaussian"))
}

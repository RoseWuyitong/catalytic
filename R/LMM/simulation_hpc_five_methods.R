# write simulatino to compare the bayesian and lmmsyn lmmridge, tau is selected by oracle prediction error
# compare, MSE of beta, prediction error, coverage rate for bayesian

# repeat 500 times for each p, and we reduce the number of p we want to see, previous p is vary from seq(10,m*n_i-5,by=5)
# now we can just use seq(10,m*n_i-10,by=10) to reduce computational burden

#  we have two bayesian methods, three point estimator MLE, ridge, cat

library(invgamma)
library(emulator)
library(MASS)
library(parallel)
numCores <- detectCores()
library(mvtnorm)
source("untility_function.R")
source("LMM_own_function.R")
source("Bayesian_LMM_function_single_sigma_u.R")


# function: return posterior mean and width of confidence interval,  coverage vector-- length is length(beta_0) contains 1 and 0
two_bayesian_given_tau <- function(syndata, mydata, tau, beta_0, ZTZ_list, ZTZ_syn_list) {
  beta_store_power_prior <- matrix(0, nc = p)
  sigma_square_store_power_prior <- c()
  sigma_u_square_store_power_prior <- c()

  beta_store_version3_same_u <- matrix(0, nc = p)
  sigma_square_store_version3_same_u <- c()
  sigma_u_square_store_version3_same_u <- c()

  beta_power_prior <- rnorm(p)
  sigma_u_square_power_prior <- 1
  sigma_square_power_prior <- 1
  random_effect_list_power_prior <- lapply(1:m, function(index) {
    rnorm(q, mean = 0, sd = sqrt(sigma_u_square_power_prior))
  })
  random_effect_syn_list_power_prior <- lapply(1:m, function(index) {
    rnorm(q, mean = 0, sd = sqrt(sigma_u_square_power_prior))
  })

  beta_version3_same_u <- beta_power_prior
  sigma_u_square_version3_same_u <- 1
  sigma_square_version3_same_u <- 1
  random_effect_list_version3_same_u <- lapply(1:m, function(index) {
    rnorm(q, mean = 0, sd = sqrt(sigma_u_square_version3_same_u))
  })


  N_MCMC <- 2000
  for (iter in 1:N_MCMC) {
    ########### power prior
    sigma_square_store_power_prior <- c(sigma_square_store_power_prior, sigma_square_power_prior)
    sigma_u_square_store_power_prior <- c(sigma_u_square_store_power_prior, sigma_u_square_power_prior)
    beta_store_power_prior <- rbind(beta_store_power_prior, beta_power_prior)
    # update beta
    beta_power_prior <- update_beta_power_prior(
      sigma_square_power_prior, sigma_u_square_power_prior, random_effect_list_power_prior,
      random_effect_syn_list_power_prior, mydata, syndata, tau
    )

    # update sigma_square
    sigma_square_power_prior <- update_sigma_square_power_prior(beta_power_prior, random_effect_list_power_prior,
      random_effect_syn_list_power_prior, mydata, syndata, tau,
      prior_rate = 0.0
    )

    # update alpha_sigma_square
    sigma_u_square_power_prior <- update_sigma_u_square_power_prior(q, random_effect_list_power_prior,
      random_effect_syn_list_power_prior, tau, M,
      prior_rate = 0.0
    )

    # update random effect list
    random_effect_list_power_prior <- update_random_effect_power_prior(beta_power_prior, q, sigma_u_square_power_prior, mydata, ZTZ_list, sigma_square_power_prior)

    # update random effect syn list
    random_effect_syn_list_power_prior <- update_random_effect_syn_power_prior(beta_power_prior, q, sigma_u_square_power_prior, syndata, ZTZ_syn_list, sigma_square_power_prior, tau, M)


    ########### version3_same_u
    sigma_square_store_version3_same_u <- c(sigma_square_store_version3_same_u, sigma_square_version3_same_u)
    sigma_u_square_store_version3_same_u <- c(sigma_u_square_store_version3_same_u, sigma_u_square_version3_same_u)
    beta_store_version3_same_u <- rbind(beta_store_version3_same_u, beta_version3_same_u)
    # update beta
    beta_version3_same_u <- update_beta_version3_same_u(
      sigma_square_version3_same_u, sigma_u_square_version3_same_u, random_effect_list_version3_same_u,
      mydata, syndata, tau
    )

    # update sigma_square
    sigma_square_version3_same_u <- update_sigma_square_version3_same_u(beta_version3_same_u, random_effect_list_version3_same_u,
      mydata, syndata, tau,
      prior_rate = 0.0
    )

    # update alpha_sigma_square
    sigma_u_square_version3_same_u <- update_sigma_u_square_version3_same_u(q, random_effect_list_version3_same_u, tau, M, prior_rate = 0.0)

    # update random effect list
    random_effect_list_version3_same_u <- update_random_effect_version3_same_u(
      beta_version3_same_u, q, sigma_u_square_version3_same_u,
      mydata, syndata, ZTZ_list, ZTZ_syn_list, sigma_square_version3_same_u, tau, M
    )
  }


  post_median_beta_power_prior <- apply(beta_store_power_prior[-(1:100), ], 2, median)
  coverage_vector_power_prior <- sapply(1:length(beta_0), function(index) {
    lower_CI <- quantile(beta_store_power_prior[-(1:100), index], 0.025)
    upper_CI <- quantile(beta_store_power_prior[-(1:100), index], 1 - 0.025)
    1 * ((beta_0[index] < upper_CI) & (beta_0[index] > lower_CI))
  })
  width_vector_power_prior <- sapply(1:length(beta_0), function(index) {
    lower_CI <- quantile(beta_store_power_prior[-(1:100), index], 0.025)
    upper_CI <- quantile(beta_store_power_prior[-(1:100), index], 1 - 0.025)
    upper_CI - lower_CI
  })

  post_median_beta_version3_same_u <- apply(beta_store_version3_same_u[-(1:100), ], 2, median)
  coverage_vector_version3_same_u <- sapply(1:length(beta_0), function(index) {
    lower_CI <- quantile(beta_store_version3_same_u[-(1:100), index], 0.025)
    upper_CI <- quantile(beta_store_version3_same_u[-(1:100), index], 1 - 0.025)
    1 * ((beta_0[index] < upper_CI) & (beta_0[index] > lower_CI))
  })
  width_vector_version3_same_u <- sapply(1:length(beta_0), function(index) {
    lower_CI <- quantile(beta_store_version3_same_u[-(1:100), index], 0.025)
    upper_CI <- quantile(beta_store_version3_same_u[-(1:100), index], 1 - 0.025)
    upper_CI - lower_CI
  })


  return(list(
    beta_power_prior = post_median_beta_power_prior, coverage_power_prior = coverage_vector_power_prior, width_power_prior = width_vector_power_prior,
    beta_version3_same_u = post_median_beta_version3_same_u, coverage_version3_same_u = coverage_vector_version3_same_u, width_version3_same_u = width_vector_version3_same_u
  ))
}







N_rep <- 250
m <- 10
n_i <- 10
p_seq <- seq(10, m * n_i - 10, by = 10)

for (rep_index in 2:(N_rep + 1)) {
  five_oracle_MSE <- matrix(0, nc = 5, nr = length(p_seq))
  five_oracle_pred <- matrix(0, nc = 5, nr = length(p_seq))
  coverage_infor_matrix_same_u <- lapply(1:length(p_seq), function(ewe) {
    rep(0, p_seq[ewe])
  })
  coverage_infor_matrix_power_prior <- lapply(1:length(p_seq), function(ewe) {
    rep(0, p_seq[ewe])
  })
  width_infor_same_u <- lapply(1:length(p_seq), function(ewe) {
    rep(0, p_seq[ewe])
  })
  width_infor_power_prior <- lapply(1:length(p_seq), function(ewe) {
    rep(0, p_seq[ewe])
  })
  for (qqq in 1:length(p_seq)) {
    p <- p_seq[qqq]
    set.seed(rep_index)
    q <- 5
    sigma_u_square <- 1
    sigma_square <- 1
    design_type <- 2
    standardize <- FALSE

    beta_0 <- generate_beta_0(p, 2, 0.5)
    Var_G <- sigma_u_square * diag(q)


    X_list <- lapply(1:m, function(index) {
      if (design_type == 1) {
        x <- cbind(1, generate_one_group_X1(n_i, p - 1, q))
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      } else {
        x <- cbind(1, generate_one_group_X2(n_i, p - 1, q))
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      }
    })

    n_i_test <- 30

    X_list_test <- lapply(1:m, function(index) {
      if (design_type == 1) {
        x <- cbind(1, generate_one_group_X1(n_i_test, p - 1, q))
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      } else {
        x <- cbind(1, generate_one_group_X2(n_i_test, p - 1, q))
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      }
    })

    Z_list <- lapply(1:m, function(index) {
      if (design_type == 1) {
        x <- generate_one_group_Z1(n_i, p, q)
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      } else {
        x <- generate_one_group_Z2(n_i, p, q)
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      }
    })

    Z_list_test <- lapply(1:m, function(index) {
      if (design_type == 1) {
        x <- generate_one_group_Z1(n_i_test, p, q)
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      } else {
        x <- generate_one_group_Z2(n_i_test, p, q)
        if (standardize) {
          meanx <- apply(x[, -1], 2, mean)
          sdx <- apply(x[, -1], 2, sd)
          x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
        }
        return(x)
      }
    })

    random_effect <- mvrnorm(m, rep(0, q), Var_G)

    Y_list <- lapply(1:m, function(index) {
      X_list[[index]] %*% beta_0 + Z_list[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i), sigma_square * diag(n_i))
    })

    Y_list_test <- lapply(1:m, function(index) {
      X_list_test[[index]] %*% beta_0 + Z_list_test[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i_test), sigma_square * diag(n_i_test))
    })

    mydata <- list(X_list = X_list, Z_list = Z_list, Y_list = Y_list)


    Y <- rep(0, m * n_i)
    for (i in 1:m) {
      Y[((i - 1) * n_i + 1):(i * n_i)] <- Y_list[[i]]
    }


    Y_test <- rep(0, m * n_i_test)
    for (i in 1:m) {
      Y_test[((i - 1) * n_i_test + 1):(i * n_i_test)] <- Y_list_test[[i]]
    }

    X <- do.call(rbind, X_list)
    X_test <- do.call(rbind, X_list_test)

    Z <- matrix(0, nr = m * n_i, nc = m * q)
    for (i in 1:m) {
      Z[((1:n_i) + n_i * (i - 1)), ((1:q) + q * (i - 1))] <- Z_list[[i]]
    }

    M <- 40

    X_list_syn <- lapply(1:m, function(index) {
      IRS_M_each_people_X(X_list[[index]], M)
    })
    Z_list_syn <- lapply(1:m, function(index) {
      IRS_M_each_people_Z(Z_list[[index]], M)
    })

    X_syn <- do.call(rbind, X_list_syn)
    Z_syn <- matrix(0, nr = m * M, nc = m * q)
    for (i in 1:m) {
      Z_syn[((1:M) + M * (i - 1)), ((1:q) + q * (i - 1))] <- Z_list_syn[[i]]
    }

    ZTZ_list <- lapply(1:m, function(i) {
      crossprod(Z_list[[i]])
    })
    ZTZ_syn_list <- lapply(1:m, function(i) {
      crossprod(Z_list_syn[[i]])
    })

    intercept_only_model <- LMM_MLE_Oct4(matrix(X[, 1], nc = 1), Y, Z)
    var_sub <- intercept_only_model$var_par
    Y_syn <- rmvnorm(1, rep(intercept_only_model$beta, M * m), sigma = var_sub[1] * diag(M * m) + var_sub[2] * tcrossprod(Z_syn))
    Y_syn <- as.numeric(Y_syn)

    grp <- rep(1:m, each = M)
    Y_list_syn <- split(Y_syn, grp)

    syndata <- list(X_list = X_list_syn, Z_list = Z_list_syn, Y_list = Y_list_syn)

    NUM_TUNING <- 50
    tau_seq <- 10^(seq(-2, 2, length.out = NUM_TUNING))
    lambda_seq_ridge <- 10^(seq(-3, 2, length.out = NUM_TUNING))
    rep_lmmridge_syn <- function(index) {
      # print(index)
      ownRidge <- LMM_ridge_Oct4(X, Y, Z, lambda = lambda_seq_ridge[index])
      ownCat <- LMM_cat_Oct4(X, Y, Z, X_syn, Y_syn, Z_syn, tau = tau_seq[index])

      matrix(c(
        norm(beta_0 - ownRidge$beta, "2")^2, lambda_seq_ridge[index], norm(beta_0 - ownCat$beta, "2")^2, tau_seq[index],
        norm(Y_test - X_test %*% (ownRidge$beta), "2")^2, lambda_seq_ridge[index], norm(Y_test - X_test %*% (ownCat$beta), "2")^2, tau_seq[index]
      ), nr = 2, byrow = T)
    }

    error_list_LMM <- mclapply(1:NUM_TUNING, rep_lmmridge_syn, mc.cores = numCores)

    error_frame_LMM_MSE <- data.frame()
    error_frame_LMM_pred <- data.frame()

    for (ii in 1:NUM_TUNING) {
      error_frame_LMM_MSE <- rbind(error_frame_LMM_MSE, error_list_LMM[[ii]][1, ])
      error_frame_LMM_pred <- rbind(error_frame_LMM_pred, error_list_LMM[[ii]][2, ])
    }

    # find the oracle tau based on the prediction error
    tau_seq <- error_frame_LMM_pred[, 4]
    oracle_tau_for_bayes <- tau_seq[which.min(error_frame_LMM_pred[, 3])]



    # write a function to take beta_0 tau mydata and syndata to produce posterior mean and coverage information

    two_baye_result <- two_bayesian_given_tau(syndata, mydata, tau = oracle_tau_for_bayes, beta_0, ZTZ_list, ZTZ_syn_list)

    # list(beta_power_prior=post_median_beta_power_prior,coverage_power_prior=coverage_vector_power_prior,width_power_prior=width_vector_power_prior,
    #     beta_version3_same_u=post_median_beta_version3_same_u,coverage_version3_same_u=coverage_vector_version3_same_u,width_version3_same_u=width_vector_version3_same_u))
    # coverage_infor_matrix_same_u=lapply(1:length(p_seq),function(ewe){matrix(0,nc=p_seq[ewe],nr=length(p_seq))})
    coverage_infor_matrix_same_u[[qqq]] <- two_baye_result$coverage_version3_same_u
    coverage_infor_matrix_power_prior[[qqq]] <- two_baye_result$coverage_power_prior

    # width_infor_same_u=lapply(1:length(p_seq),function(ewe){rep(0,p_seq[ewe])})
    # width_infor_power_prior=lapply(1:length(p_seq),function(ewe){rep(0,p_seq[ewe])})
    width_infor_same_u[[qqq]] <- two_baye_result$width_version3_same_u
    width_infor_power_prior[[qqq]] <- two_baye_result$width_power_prior

    ownMLE <- LMM_MLE_Oct4(X, Y, Z)
    five_oracle_MSE[qqq, ] <- c(norm(beta_0 - ownMLE$beta, "2")^2, min(error_frame_LMM_MSE[, 1]), min(error_frame_LMM_MSE[, 3]), norm(beta_0 - two_baye_result$beta_power_prior, "2")^2, norm(beta_0 - two_baye_result$beta_version3_same_u, "2")^2)
    five_oracle_pred[qqq, ] <- c(norm(Y_test - X_test %*% (ownMLE$beta), "2")^2, min(error_frame_LMM_pred[, 1]), min(error_frame_LMM_pred[, 3]), norm(Y_test - X_test %*% (two_baye_result$beta_power_prior), "2")^2, norm(Y_test - X_test %*% (two_baye_result$beta_version3_same_u), "2")^2)
  }

  five_oracle_MSE <- as.data.frame(five_oracle_MSE)
  five_oracle_pred <- as.data.frame(five_oracle_pred)
  colnames(five_oracle_MSE) <- c("MLE", "ridge", "Cat", "power_prior", "ver3_same_u")
  colnames(five_oracle_pred) <- c("MLE", "ridge", "Cat", "power_prior", "ver3_same_u")
  five_oracle_pred$p <- p_seq
  five_oracle_MSE$p <- p_seq
  save(five_oracle_MSE, file = paste0("five_oracle_MSE", rep_index))
  save(five_oracle_pred, file = paste0("five_oracle_pred", rep_index))
  save(width_infor_same_u, file = paste0("width_same_u", rep_index))
  save(width_infor_power_prior, file = paste0("width_power_prior", rep_index))
  save(coverage_infor_matrix_same_u, file = paste0("coverage_same_u", rep_index))
  save(coverage_infor_matrix_power_prior, file = paste0("coverage_power_prior", rep_index))
}

# function used in Feb 14 script





# function: compute MLE and estimation of sigma by gradient descent
# input: dataset X_list,Z_list,Y_list
# output: beta_mle, and variance for random effect and noise
MLE_unknown_G <- function(X_list, Z_list, Y_list) {
  p <- ncol(X_list[[1]])
  MLE_result <- GD_MLE(X_list, Y_list, Z_list, rep(0, p + 2), 0.01)
  return(list(betahat = MLE_result$beta, sigma_square = MLE_result$sigma_sq, sigma_u_square = MLE_result$sigma_u_sq))
}
# function: return 5 datasets based on original data, each data set contain one training set and one validation set
five_fold_CV_separation <- function(X_list, Z_list, Y_list) {
  n_i <- nrow(X_list[[1]])
  m <- length(X_list)
  ind <- matrix(1:n_i, byrow = T, nr = 5)
  whole_set <- list()
  for (j in 1:5) {
    valid_X <- lapply(1:m, function(i) {
      X_i <- X_list[[i]]
      X_i[ind[j, ], ]
    })
    X_list_train <- lapply(1:m, function(i) {
      X_i <- X_list[[i]]
      X_i[-ind[j, ], ]
    })
    Z_list_train <- lapply(1:m, function(i) {
      Z_i <- Z_list[[i]]
      Z_i[-ind[j, ], ]
    })
    valid_Z <- lapply(1:m, function(i) {
      Z_i <- Z_list[[i]]
      Z_i[ind[j, ], ]
    })
    Y_list_train <- lapply(1:m, function(i) {
      Y_i <- Y_list[[i]]
      Y_i[-ind[j, ]]
    })
    valid_Y <- lapply(1:m, function(i) {
      Y_i <- Y_list[[i]]
      Y_i[ind[j, ]]
    })
    whole_set[[j]] <- list(
      valid_X = valid_X, valid_Z = valid_Z, valid_Y = valid_Y,
      X_list_train = X_list_train, Z_list_train = Z_list_train, Y_list_train = Y_list_train
    )
  }
  return(whole_set)
}

# function: cross validation for ridge estimator
# input: X,Y,Z_list, and data used for cross-validation
# output: four estimator, predict error |a-b|: no estimation of random effect, has random effect
#                         predict error (a-b)^2: no estimation of random effect, has random effect
#                       result in four tuning parameter
# only need to export estimator and estimation of variance, the prediction is done in another big function, since all estimation need to be compared on same testing set
library(parallel)
n.cores <- detectCores()
five_fold_CV_ridge_unknown_G <- function(X_list, Z_list, Y_list, whole_set) {
  lambda_seq <- 10^seq(-6, 2, by = .35)
  p <- ncol(X_list[[1]])
  rep_f <- function(lambda_) {
    m <- length(whole_set[[1]]$X_list_train)
    library(MASS)
    library(psych)
    GD_ridge <- function(X_list, Y_list, Z_list, x_init, lambda, step_size) {
      tol <- 10^(-4)
      momentum <- 0.5
      p <- length(x_init) - 2
      gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
        m <- length(X_list)
        n_i <- nrow(X_list[[1]])
        gd <- apply(sapply(1:m, function(i) {
          Xi <- X_list[[i]]
          Yi <- Y_list[[i]]
          Zi <- Z_list[[i]]
          Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))
          component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)

          component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
          component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
          return(c(component1, component2, component3))
        }), 1, sum)
        return(-gd + lambda * c(beta, 0, 0))
      }
      delta_x_old <- 0
      x_old <- x_init
      delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
      x_new <- x_old - delta_x_new
      cnt <- 1
      iter <- 1
      while (norm(x_new - x_old, "2") > tol & iter < 50000) {
        delta_x_old <- delta_x_new
        x_old <- x_new
        delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
        x_new <- x_old - delta_x_new
        cnt <- cnt + 1
        iter <- iter + 1
        if (cnt > 500) {
          step_size <- 0.5 * step_size
          cnt <- cnt - 500
        }
        # cat("iteration: ",iter,"\n")
      }
      # cat("iteration: ",iter,"\n")
      return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2]), iteration = iter))
    }
    # in this rep_f ,we return the error for four kind of estimation method we mention above, that is to say, we need to use parLapply
    # each parallel list contain four CV error at current lambda_
    # Luckily, we only to use GD_ridge five times in this parallel thread
    error_metric <- rep(0, 5) # this error metric will be average of five validation error
    #                     abs_no_random_effect, L2_no_random_effect, abs_has_random effect, L2_has_random_effect
    for (j in 1:5) {
      current_dat <- whole_set[[j]]
      X_list_cur <- current_dat$X_list_train
      Z_list_cur <- current_dat$Z_list_train
      Y_list_cur <- current_dat$Y_list_train
      p <- ncol(X_list_cur[[1]])
      m <- length(X_list_cur)
      n_i <- nrow(X_list_cur[[1]])
      q <- ncol(Z_list_cur[[1]])
      cur_GD_ridge <- GD_ridge(X_list_cur, Y_list_cur, Z_list_cur, rep(0, p + 2), lambda = lambda_, 0.01)

      cur_beta_j <- cur_GD_ridge$beta
      # cat(cur_beta_j,"\n")
      Var_G_cur <- diag(q) * cur_GD_ridge$sigma_u_sq
      V_list_cur <- lapply(1:m, function(index) {
        cur_GD_ridge$sigma_u_sq * diag(n_i) + Z_list_cur[[index]] %*% Var_G_cur %*% t(Z_list_cur[[index]])
      })
      mixed_effect_pred_BLUE <- lapply(1:m, function(index) {
        Var_G_cur %*% t(Z_list_cur[[index]]) %*% ginv(V_list_cur[[index]]) %*% (Y_list_cur[[index]] - X_list_cur[[index]] %*% cur_beta_j)
      })
      prediction_function_has_random_effect_inside_CV_L2 <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat + Z_list_test[[index]] %*% mixed_effect_pred_BLUE[[index]]
          norm(pp - Y_list_test[[index]], "2")^2
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_has_random_effect_inside_CV_abs <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat + Z_list_test[[index]] %*% mixed_effect_pred_BLUE[[index]]
          sum(abs(pp - Y_list_test[[index]]))
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_no_random_effect_inside_CV_L2 <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat
          norm(pp - Y_list_test[[index]], "2")^2
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_no_random_effect_inside_CV_abs <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat
          sum(abs(pp - Y_list_test[[index]]))
        })
        return(sum(Y_list_pred_error) / (m))
      }
      error_metric <- error_metric + c(
        prediction_function_no_random_effect_inside_CV_abs(cur_beta_j),
        prediction_function_no_random_effect_inside_CV_L2(cur_beta_j),
        prediction_function_has_random_effect_inside_CV_abs(cur_beta_j),
        prediction_function_has_random_effect_inside_CV_L2(cur_beta_j), cur_GD_ridge$iteration
      )
    }
    error_metric <- error_metric / 5
    return(c(error_metric, lambda_))
  }

  clust <- makeCluster(n.cores)
  clusterExport(clust, "whole_set", envir = environment())
  result_all_lambda <- parLapply(clust, lambda_seq, rep_f)
  stopCluster(clust)
  #                     abs_no_random_effect, L2_no_random_effect, abs_has_random effect, L2_has_random_effect
  matrix_version_result_all_lambda <- matrix(0, nr = length(lambda_seq), nc = 6)
  for (i in 1:length(lambda_seq)) {
    matrix_version_result_all_lambda[i, ] <- result_all_lambda[[i]]
  }
  cat("iteration for each lambda is: ", matrix_version_result_all_lambda[, 5], "\n")
  best_lambda_abs_no_random_effect <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 1]), 6]
  best_lambda_L2_no_random_effect <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 2]), 6]
  best_lambda_abs_has_random_effect <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 3]), 6]
  best_lambda_L2_has_random_effect <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 4]), 6]
  result_ridge_abs_no_random_effect <- GD_ridge(X_list, Y_list, Z_list, rep(0, p + 2), best_lambda_abs_no_random_effect, 0.01)
  result_ridge_L2_no_random_effect <- GD_ridge(X_list, Y_list, Z_list, rep(0, p + 2), best_lambda_L2_no_random_effect, 0.01)
  result_ridge_abs_has_random_effect <- GD_ridge(X_list, Y_list, Z_list, rep(0, p + 2), best_lambda_abs_has_random_effect, 0.01)
  result_ridge_L2_has_random_effect <- GD_ridge(X_list, Y_list, Z_list, rep(0, p + 2), best_lambda_L2_has_random_effect, 0.01)

  beta_ridge_abs_no_random_effect <- result_ridge_abs_no_random_effect$beta
  beta_ridge_L2_no_random_effect <- result_ridge_L2_no_random_effect$beta
  beta_ridge_abs_has_random_effect <- result_ridge_abs_has_random_effect$beta
  beta_ridge_L2_has_random_effect <- result_ridge_L2_has_random_effect$beta

  ridge_abs_no_random_effect <- list(
    betahat = beta_ridge_abs_no_random_effect, sigma_square = result_ridge_abs_no_random_effect$sigma_sq,
    sigma_u_square = result_ridge_abs_no_random_effect$sigma_u_sq, best_lambda = best_lambda_abs_no_random_effect
  )
  ridge_L2_no_random_effect <- list(
    betahat = beta_ridge_L2_no_random_effect, sigma_square = result_ridge_L2_no_random_effect$sigma_sq,
    sigma_u_square = result_ridge_L2_no_random_effect$sigma_u_sq, best_lambda = best_lambda_L2_no_random_effect
  )
  ridge_abs_has_random_effect <- list(
    betahat = beta_ridge_abs_has_random_effect, sigma_square = result_ridge_abs_has_random_effect$sigma_sq,
    sigma_u_square = result_ridge_abs_has_random_effect$sigma_u_sq, best_lambda = best_lambda_abs_has_random_effect
  )
  ridge_L2_has_random_effect <- list(
    betahat = beta_ridge_L2_has_random_effect, sigma_square = result_ridge_L2_has_random_effect$sigma_sq,
    sigma_u_square = result_ridge_L2_has_random_effect$sigma_u_sq, best_lambda = best_lambda_L2_has_random_effect
  )
  return(list(
    ridge_abs_no_random_effect = ridge_abs_no_random_effect, ridge_L2_no_random_effect = ridge_L2_no_random_effect,
    ridge_abs_has_random_effect = ridge_abs_has_random_effect, ridge_L2_has_random_effect = ridge_L2_has_random_effect
  ))
}


# function: cross validation for catalytic estimator
# input: X,Y,Z_list and data used for cross validation
# output: four estimator, predict error |a-b|: no estimation of random effect, has random effect
#                         predict error (a-b)^2: no estimation of random effect, has random effect
#                       result in four tuning parameter
# create synthetic X and Z first, synthetic Y will depends on you initial beta_0_hat, e.g submodel is first column and last column(abandon)
# !!!!!!! attention only synthetic response is mean is consider for simplicity
five_fold_CV_cat_unknown_G <- function(X_list, Z_list, Y_list, whole_set, M.syn) {
  tau_seq <- 10^seq(-5, 3, by = .35)
  tau_seq <- tau_seq[tau_seq < M.syn]
  M <- M.syn
  p <- ncol(X_list[[1]])
  m <- length(X_list)
  X_list_syn <- lapply(1:m, function(index) {
    IRS_M_each_people_X(X_list[[index]], M)
  })
  Z_list_syn <- lapply(1:m, function(index) {
    IRS_M_each_people_Z(Z_list[[index]], M)
  })
  Y_list_syn <- lapply(1:m, function(index) {
    rep(mean(Y_list[[index]]), M)
  })
  rep_f2 <- function(tau_) {
    m <- length(whole_set[[1]]$X_list_train)
    library(MASS)
    library(psych)
    GD_cat <- function(X_list, Y_list, Z_list, x_init, tau_, X_list_syn, Y_list_syn, Z_list_syn, step_size) {
      tol <- 10^(-4)
      p <- length(x_init) - 2
      m <- length(X_list)
      gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
        m <- length(X_list)
        n_i <- nrow(X_list[[1]])
        gd <- apply(sapply(1:m, function(i) {
          Xi <- X_list[[i]]
          Yi <- Y_list[[i]]
          Zi <- Z_list[[i]]
          Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

          component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
          component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
          component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
          return(c(component1, component2, component3))
        }), 1, sum)
        return(-gd)
      }
      gradient_cat <- function(beta, A, B, X_list, Y_list, Z_list) {
        m <- length(X_list)
        n_i <- nrow(X_list[[1]])
        gd <- apply(sapply(1:m, function(i) {
          Xi <- X_list[[i]]
          Yi <- Y_list[[i]]
          Zi <- Z_list[[i]]
          Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

          component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)

          return(c(component1, 0, 0))
        }), 1, sum)
        return(-gd)
      }
      x_old <- x_init
      M.syn <- nrow(X_list_syn[[1]])
      gradient_merge <- gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + tau_ / M.syn * gradient_cat(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list_syn, Y_list_syn, Z_list_syn)
      x_new <- x_old - step_size * gradient_merge
      cnt <- 1
      iter <- 1
      while (norm(x_new - x_old, "2") > tol & iter < 5000) {
        x_old <- x_new
        gradient_merge <- gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + tau_ / M.syn * gradient_cat(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list_syn, Y_list_syn, Z_list_syn)
        x_new <- x_old - step_size * gradient_merge
        cnt <- cnt + 1
        iter <- iter + 1
        if (cnt > 200) {
          step_size <- 0.8 * step_size
          cnt <- cnt - 200
        }
        # cat("iteration: ",iter,"\n")
      }
      return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2]), iteratioin = iter))
    }

    # in this rep_f ,we return the error for four kind of estimation method we mention above, that is to say, we need to use parLapply
    # each parallel list contain four CV error at current lambda_
    # Luckily, we only to use GD_cat five times in this parallel thread
    # noticed that our synthetic response is just the mean, we need to change to another Y_list_syn
    error_metric <- rep(0, 5) # this error metric will be average of five validation error
    #                     abs_no_random_effect, L2_no_random_effect, abs_has_random effect, L2_has_random_effect
    for (j in 1:5) {
      current_dat <- whole_set[[j]]
      X_list_cur <- current_dat$X_list_train
      Z_list_cur <- current_dat$Z_list_train
      Y_list_cur <- current_dat$Y_list_train
      p <- ncol(X_list_cur[[1]])
      m <- length(X_list_cur)
      n_i <- nrow(X_list_cur[[1]])
      q <- ncol(Z_list_cur[[1]])
      cur_GD_cat <- GD_cat(X_list_cur, Y_list_cur, Z_list_cur, rep(0, p + 2), tau_, X_list_syn, Y_list_syn, Z_list_syn, 0.01)

      cur_beta_j <- cur_GD_cat$beta
      # cat(cur_beta_j,"\n")
      Var_G_cur <- diag(q) * cur_GD_cat$sigma_u_sq
      V_list_cur <- lapply(1:m, function(index) {
        cur_GD_cat$sigma_u_sq * diag(n_i) + Z_list_cur[[index]] %*% Var_G_cur %*% t(Z_list_cur[[index]])
      })
      mixed_effect_pred_BLUE <- lapply(1:m, function(index) {
        Var_G_cur %*% t(Z_list_cur[[index]]) %*% ginv(V_list_cur[[index]]) %*% (Y_list_cur[[index]] - X_list_cur[[index]] %*% cur_beta_j)
      })
      prediction_function_has_random_effect_inside_CV_L2 <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat + Z_list_test[[index]] %*% mixed_effect_pred_BLUE[[index]]
          norm(pp - Y_list_test[[index]], "2")^2
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_has_random_effect_inside_CV_abs <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat + Z_list_test[[index]] %*% mixed_effect_pred_BLUE[[index]]
          sum(abs(pp - Y_list_test[[index]]))
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_no_random_effect_inside_CV_L2 <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat
          norm(pp - Y_list_test[[index]], "2")^2
        })
        return(sum(Y_list_pred_error) / (m))
      }
      prediction_function_no_random_effect_inside_CV_abs <- function(beta_hat) {
        X_list_test <- current_dat$valid_X
        Z_list_test <- current_dat$valid_Z
        Y_list_test <- current_dat$valid_Y
        Y_list_pred_error <- sapply(1:m, function(index) {
          pp <- X_list_test[[index]] %*% beta_hat
          sum(abs(pp - Y_list_test[[index]]))
        })
        return(sum(Y_list_pred_error) / (m))
      }
      error_metric <- error_metric + c(
        prediction_function_no_random_effect_inside_CV_abs(cur_beta_j),
        prediction_function_no_random_effect_inside_CV_L2(cur_beta_j),
        prediction_function_has_random_effect_inside_CV_abs(cur_beta_j),
        prediction_function_has_random_effect_inside_CV_L2(cur_beta_j), cur_GD_cat$iteration
      )
    }
    error_metric <- error_metric / 5
    return(c(error_metric, tau_))
  }
  clust <- makeCluster(n.cores)
  clusterExport(clust, varlist = c("whole_set", "X_list_syn", "Z_list_syn", "Y_list_syn", "M.syn"), envir = environment())
  result_all_tau <- parLapply(clust, tau_seq, rep_f2)
  stopCluster(clust)
  #                     abs_no_random_effect, L2_no_random_effect, abs_has_random effect, L2_has_random_effect
  matrix_version_result_all_tau <- matrix(0, nr = length(tau_seq), nc = 6)
  # cat(length(result_all_tau[[1]]),dim(matrix_version_result_all_tau),length(tau_seq),"\n")
  for (i in 1:length(tau_seq)) {
    matrix_version_result_all_tau[i, ] <- result_all_tau[[i]]
  }
  cat("iteration for each tau is: ", matrix_version_result_all_tau[, 5], "\n")
  best_tau_abs_no_random_effect <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 1]), 6]
  best_tau_L2_no_random_effect <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 2]), 6]
  best_tau_abs_has_random_effect <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 3]), 6]
  best_tau_L2_has_random_effect <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 4]), 6]
  result_cat_abs_no_random_effect <- GD_cat(X_list, Y_list, Z_list, rep(0, p + 2), best_tau_abs_no_random_effect, X_list_syn, Y_list_syn, Z_list_syn, 0.01)
  result_cat_L2_no_random_effect <- GD_cat(X_list, Y_list, Z_list, rep(0, p + 2), best_tau_L2_no_random_effect, X_list_syn, Y_list_syn, Z_list_syn, 0.01)
  result_cat_abs_has_random_effect <- GD_cat(X_list, Y_list, Z_list, rep(0, p + 2), best_tau_abs_has_random_effect, X_list_syn, Y_list_syn, Z_list_syn, 0.01)
  result_cat_L2_has_random_effect <- GD_cat(X_list, Y_list, Z_list, rep(0, p + 2), best_tau_L2_has_random_effect, X_list_syn, Y_list_syn, Z_list_syn, 0.01)

  beta_cat_abs_no_random_effect <- result_cat_abs_no_random_effect$beta
  beta_cat_L2_no_random_effect <- result_cat_L2_no_random_effect$beta
  beta_cat_abs_has_random_effect <- result_cat_abs_has_random_effect$beta
  beta_cat_L2_has_random_effect <- result_cat_L2_has_random_effect$beta

  cat_abs_no_random_effect <- list(
    betahat = beta_cat_abs_no_random_effect, sigma_square = result_cat_abs_no_random_effect$sigma_sq,
    sigma_u_square = result_cat_abs_no_random_effect$sigma_u_sq, best_tau = best_tau_abs_no_random_effect
  )
  cat_L2_no_random_effect <- list(
    betahat = beta_cat_L2_no_random_effect, sigma_square = result_cat_L2_no_random_effect$sigma_sq,
    sigma_u_square = result_cat_L2_no_random_effect$sigma_u_sq, best_tau = best_tau_L2_no_random_effect
  )
  cat_abs_has_random_effect <- list(
    betahat = beta_cat_abs_has_random_effect, sigma_square = result_cat_abs_has_random_effect$sigma_sq,
    sigma_u_square = result_cat_abs_has_random_effect$sigma_u_sq, best_tau = best_tau_abs_has_random_effect
  )
  cat_L2_has_random_effect <- list(
    betahat = beta_cat_L2_has_random_effect, sigma_square = result_cat_L2_has_random_effect$sigma_sq,
    sigma_u_square = result_cat_L2_has_random_effect$sigma_u_sq, best_tau = best_tau_L2_has_random_effect
  )
  return(list(
    cat_abs_no_random_effect = cat_abs_no_random_effect, cat_L2_no_random_effect = cat_L2_no_random_effect,
    cat_abs_has_random_effect = cat_abs_has_random_effect, cat_L2_has_random_effect = cat_L2_has_random_effect
  ))
}


one_time_simulation_MLE_ridge_cat_oracle_diff_init_CV <- function(m, n_i, beta_0, q, sigma_u_square, sigma_square, M.syn, design_type, rep_index) {
  p <- length(beta_0)
  Var_G <- sigma_u_square * diag(q)
  X_list <- lapply(1:m, function(index) {
    if (design_type == 1) {
      return(generate_one_group_X1(n_i, p, q))
    } else {
      return(generate_one_group_X2(n_i, p, q))
    }
  })
  Z_list <- lapply(1:m, function(index) {
    if (design_type == 1) {
      return(generate_one_group_Z1(n_i, p, q))
    } else {
      return(generate_one_group_Z2(n_i, p, q))
    }
  })
  V_list <- lapply(1:m, function(index) {
    sigma_square^2 * diag(n_i) + Z_list[[index]] %*% Var_G %*% t(Z_list[[index]])
  })
  random_effect <- mvrnorm(m, rep(0, q), Var_G)

  Y_list <- lapply(1:m, function(index) {
    X_list[[index]] %*% beta_0 + Z_list[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i), diag(n_i))
  })

  whole_set <- five_fold_CV_separation(X_list, Z_list, Y_list)
  synthetic_result <- five_fold_CV_cat_unknown_G(X_list, Z_list, Y_list, whole_set, M.syn = 50)
  cat(rep_index, " is okey for catalytic \n")
  ridge_result <- five_fold_CV_ridge_unknown_G(X_list, Z_list, Y_list, whole_set)
  cat(rep_index, " is okey for ridge \n")

  MLE_result <- MLE_unknown_G(X_list, Z_list, Y_list)
  cat(rep_index, " is okey for MLE \n")
  original_data <- list(X = X_list, Z = Z_list, Y = Y_list, true_beta = beta_0, random_effect = random_effect)
  name_ori <- paste0("original_data", rep_index)
  name_MLE <- paste0("MLE", rep_index)
  name_ridge <- paste0("ridge", rep_index)
  name_cat <- paste0("cat", rep_index)
  save(ridge_result, file = name_ridge)
  save(MLE_result, file = name_MLE)
  save(synthetic_result, file = name_cat)
  save(original_data, file = name_ori)
}

ten_time_repeat <- function(m = 10, n_i = 10, p = 50, q = 10, sigma_u_square = 1, sigma_square = 1, beta_0_prop = 0.75, M.syn = 50, design_type = 1) {
  beta_0 <- generate_beta_0(p, 2, beta_0_prop)
  for (k in 1:10) {
    one_time_simulation_MLE_ridge_cat_oracle_diff_init_CV(m, n_i, beta_0, q, sigma_u_square, sigma_square, M.syn, design_type, k)
  }
}



# function: cross validation for lmmridge, not initialized by GD
# input: X,Y,Z_list, and data used for cross-validation
# output: two estimator with their best tuning parameter ridge_abs=ridge_abs,ridge_L2=ridge_L2,MLE=beta_MLE


five_fold_CV_lmmridge <- function(X_list, Z_list, Y_list, whole_set, lambda_seq) {
  length_lambda_seq <- length(lambda_seq)
  p <- ncol(X_list[[1]])
  m <- length(X_list)
  n_i <- nrow(X_list[[1]])
  q <- ncol(Z_list[[1]])


  rep_f <- function(lambda_) {
    # m=length(whole_set[[1]]$X_list_train)
    # library(MASS)
    # library(psych)
    # in this rep_f ,we return the error for four kind of estimation method we mention above, that is to say, we need to use parLapply

    error_metric <- rep(0, 2) # this error metric will be average of five validation error
    #                      abs_has_random effect, L2_has_random_effect
    for (j in 1:5) {
      current_dat <- whole_set[[j]]
      X_list_cur <- current_dat$X_list_train
      Z_list_cur <- current_dat$Z_list_train
      Y_list_cur <- current_dat$Y_list_train
      X_list_test_cur <- current_dat$valid_X
      Y_list_test_cur <- current_dat$valid_Y
      Z_list_test_cur <- current_dat$valid_Z
      p_cur <- ncol(X_list_cur[[1]])
      m_cur <- length(X_list_cur)
      n_i_cur <- nrow(X_list_cur[[1]])
      q_cur <- ncol(Z_list_cur[[1]])
      n_i_valid <- nrow(X_list_test_cur[[1]])
      mydata <- matrix(0, nc = p_cur + q_cur + 2, nr = m_cur * n_i_cur)

      for (i in 1:m_cur) {
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), 1] <- Y_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), 2:(p_cur + 1)] <- X_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), (p_cur + 2):(p_cur + 1 + q_cur)] <- Z_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), (p_cur + q_cur + 2)] <- i
      }

      mydata <- as.data.frame(mydata)
      colnames(mydata) <- c("Y", paste0("X", 1:p_cur), paste0("Z", 1:q_cur), "group")
      mydata_valid <- matrix(0, nc = p_cur + q_cur + 2, nr = m_cur * n_i_valid)
      for (i in 1:m_cur) {
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), 1] <- Y_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), 2:(p_cur + 1)] <- X_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), (p_cur + 2):(p_cur + 1 + q_cur)] <- Z_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), (p_cur + q_cur + 2)] <- i
      }
      mydata_valid <- as.data.frame(mydata_valid)
      colnames(mydata_valid) <- c("Y", paste0("X", 1:p_cur), paste0("Z", 1:q_cur), "group")


      lmm_ridge_cur <- lmmridge.pred_testdata(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
        mydata$group,
        lambda = lambda_, startValue = 2,
        pdMat = "pdIdent", CovOpt = "nlm", standardize = F, test_x = as.matrix(mydata_valid[, 2:(p + 1)]),
        test_y = mydata_valid$Y, test_z = as.matrix(mydata_valid[, (p + 2):(p + q + 1)]), test_grp = mydata_valid$group
      )
      residual_cur <- as.numeric(lmm_ridge_cur$residuals)


      error_metric <- error_metric + c(sum(abs(residual_cur)), sum(residual_cur^2))
    }
    error_metric <- error_metric / 5
    return(c(error_metric, lambda_))
  }
  ################################################# nonparallel version
  result_all_lambda <- list()
  length(result_all_lambda) <- length_lambda_seq
  for (k in 1:length_lambda_seq) {
    result_all_lambda[[k]] <- rep_f(lambda_seq[k])
    # cat(k,"-th lambda finish")
  }


  matrix_version_result_all_lambda <- matrix(0, nr = length(lambda_seq), nc = length(result_all_lambda[[1]]))
  for (i in 1:length(lambda_seq)) {
    matrix_version_result_all_lambda[i, ] <- result_all_lambda[[i]]
  }
  best_lambda_abs_pred <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 1]), 3]
  best_lambda_L2_pred <- matrix_version_result_all_lambda[which.min(matrix_version_result_all_lambda[, 2]), 3]
  mydata <- matrix(0, nc = p + q + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- X_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- Z_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("Y", paste0("X", 1:p), paste0("Z", 1:q), "group")
  result_ridge_abs <- lmmridge(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = best_lambda_abs_pred, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F
  )
  result_ridge_L2 <- lmmridge(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = best_lambda_L2_pred, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F
  )
  result_MLE <- lmmridge(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = 0, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F
  )
  beta_ridge_abs <- result_ridge_abs$coefficients
  beta_ridge_L2 <- result_ridge_L2$coefficients
  beta_MLE <- result_MLE$coefficients
  # we do not return the variance parameter since they are not our objective
  ridge_abs <- list(betahat = beta_ridge_abs, best_lambda = best_lambda_abs_pred)
  ridge_L2 <- list(betahat = beta_ridge_L2, best_lambda = best_lambda_L2_pred)

  return(list(
    ridge_abs = ridge_abs, ridge_L2 = ridge_L2, MLE = beta_MLE, matrix_version_result_all_lambda = matrix_version_result_all_lambda,
    random_abs = result_ridge_abs$random, random_L2 = result_ridge_L2$random
  ))
}


#  function: cross validation for lmm.syn, not initialized by GD
# input: X,Y,Z_list, and data used for cross-validation ,tau_seq and syn_data
# output: two estimator with their best tuning parameter syn_abs,syn_L2,MLE=beta_MLE

five_fold_CV_lmmsyn <- function(X_list, Z_list, Y_list, whole_set, tau_seq, my_syndata) {
  length_tau_seq <- length(tau_seq)
  p <- ncol(X_list[[1]])
  m <- length(X_list)
  n_i <- nrow(X_list[[1]])
  q <- ncol(Z_list[[1]])

  rep_f <- function(tau_) {
    # m=length(whole_set[[1]]$X_list_train)
    # library(MASS)
    # library(psych)
    # in this rep_f ,we return the error for four kind of estimation method we mention above, that is to say, we need to use parLapply

    error_metric <- rep(0, 2) # this error metric will be average of five validation error
    #                      abs_has_random effect, L2_has_random_effect
    for (j in 1:5) {
      current_dat <- whole_set[[j]]
      X_list_cur <- current_dat$X_list_train
      Z_list_cur <- current_dat$Z_list_train
      Y_list_cur <- current_dat$Y_list_train
      X_list_test_cur <- current_dat$valid_X
      Y_list_test_cur <- current_dat$valid_Y
      Z_list_test_cur <- current_dat$valid_Z
      p_cur <- ncol(X_list_cur[[1]])
      m_cur <- length(X_list_cur)
      n_i_cur <- nrow(X_list_cur[[1]])
      q_cur <- ncol(Z_list_cur[[1]])
      n_i_valid <- nrow(X_list_test_cur[[1]])
      mydata <- matrix(0, nc = p_cur + q_cur + 2, nr = m_cur * n_i_cur)

      for (i in 1:m_cur) {
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), 1] <- Y_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), 2:(p_cur + 1)] <- X_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), (p_cur + 2):(p_cur + 1 + q_cur)] <- Z_list_cur[[i]]
        mydata[((i - 1) * n_i_cur + 1):(i * n_i_cur), (p_cur + q_cur + 2)] <- i
      }

      mydata <- as.data.frame(mydata)
      colnames(mydata) <- c("Y", paste0("X", 1:p_cur), paste0("Z", 1:q_cur), "group")
      mydata_valid <- matrix(0, nc = p_cur + q_cur + 2, nr = m_cur * n_i_valid)
      for (i in 1:m_cur) {
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), 1] <- Y_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), 2:(p_cur + 1)] <- X_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), (p_cur + 2):(p_cur + 1 + q_cur)] <- Z_list_test_cur[[i]]
        mydata_valid[((i - 1) * n_i_valid + 1):(i * n_i_valid), (p_cur + q_cur + 2)] <- i
      }
      mydata_valid <- as.data.frame(mydata_valid)
      colnames(mydata_valid) <- c("Y", paste0("X", 1:p_cur), paste0("Z", 1:q_cur), "group")



      lmm_syn_cur <- lmmsyn.pred_testdata(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
        mydata$group,
        lambda = 0, startValue = 2,
        pdMat = "pdIdent", CovOpt = "nlm", standardize = F, my_syndata = my_syndata, tau_syn = tau_, test_x = as.matrix(mydata_valid[, 2:(p + 1)]),
        test_y = mydata_valid$Y, test_z = as.matrix(mydata_valid[, (p + 2):(p + q + 1)]), test_grp = mydata_valid$group
      )
      residual_cur <- as.numeric(lmm_syn_cur$residuals)


      error_metric <- error_metric + c(sum(abs(residual_cur)), sum(residual_cur^2))
    }
    error_metric <- error_metric / 5
    return(c(error_metric, tau_))
  }

  ################################################# nonparallel version
  result_all_tau <- list()
  length(result_all_tau) <- length_tau_seq
  for (k in 1:length_tau_seq) {
    result_all_tau[[k]] <- rep_f(tau_seq[k])
  }
  matrix_version_result_all_tau <- matrix(0, nr = length(tau_seq), nc = length(result_all_tau[[1]]))
  for (i in 1:length(tau_seq)) {
    matrix_version_result_all_tau[i, ] <- result_all_tau[[i]]
  }
  best_tau_abs_pred <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 1]), 3]
  best_tau_L2_pred <- matrix_version_result_all_tau[which.min(matrix_version_result_all_tau[, 2]), 3]
  mydata <- matrix(0, nc = p + q + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- X_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- Z_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("Y", paste0("X", 1:p), paste0("Z", 1:q), "group")

  result_syn_abs <- lmmsyn(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = 0, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F, my_syndata = my_syndata, tau_syn = best_tau_abs_pred
  )

  result_syn_L2 <- lmmsyn(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = 0, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F, my_syndata = my_syndata, tau_syn = best_tau_L2_pred
  )

  result_MLE <- lmmsyn(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = 0, startValue = 2,
    pdMat = "pdIdent", CovOpt = "nlm", standardize = F, my_syndata = my_syndata, tau_syn = 0
  )

  beta_syn_abs <- result_syn_abs$coefficients
  beta_syn_L2 <- result_syn_L2$coefficients
  beta_MLE <- result_MLE$coefficients
  # we do not return the variance parameter since they are not our objective
  syn_abs <- list(betahat = beta_syn_abs, best_tau = best_tau_abs_pred)
  syn_L2 <- list(betahat = beta_syn_L2, best_tau = best_tau_L2_pred)

  return(list(
    syn_abs = syn_abs, syn_L2 = syn_L2, MLE = beta_MLE, matrix_version_result_all_tau = matrix_version_result_all_tau,
    random_abs = result_syn_abs$random, random_L2 = result_syn_L2$random
  ))
}






































GD_ridge <- function(X_list, Y_list, Z_list, x_init, lambda, step_size) {
  tol <- 10^(-4)
  momentum <- 0.5
  p <- length(x_init) - 2
  gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd + lambda * c(beta, 0, 0))
  }
  delta_x_old <- 0
  x_old <- x_init
  delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
  x_new <- x_old - delta_x_new
  cnt <- 1
  iter <- 1
  while (norm(x_new - x_old, "2") > tol & iter < 50000) {
    delta_x_old <- delta_x_new
    x_old <- x_new
    delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
    x_new <- x_old - delta_x_new
    cnt <- cnt + 1
    iter <- iter + 1
    if (cnt > 200) {
      step_size <- 0.8 * step_size
      cnt <- cnt - 200
    }
    #  cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2])))
}
GD_cat <- function(X_list, Y_list, Z_list, x_init, tau, X_list_syn, Y_list_syn, Z_list_syn, step_size) {
  tol <- 10^(-4)
  p <- length(x_init) - 2
  gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd)
  }
  gradient_cat <- function(beta, A, B, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)

      return(c(component1, 0, 0))
    }), 1, sum)
    return(-gd)
  }
  x_old <- x_init
  M.syn <- nrow(X_list_syn[[1]])
  gradient_merge <- gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + tau / M.syn * gradient_cat(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list_syn, Y_list_syn, Z_list_syn)
  x_new <- x_old - step_size * gradient_merge
  cnt <- 1
  iter <- 1
  while (norm(x_new - x_old, "2") > tol & iter < 5000) {
    x_old <- x_new
    gradient_merge <- gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + tau / M.syn * gradient_cat(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list_syn, Y_list_syn, Z_list_syn)
    x_new <- x_old - step_size * gradient_merge
    cnt <- cnt + 1
    iter <- iter + 1
    if (cnt > 200) {
      step_size <- 0.8 * step_size
      cnt <- cnt - 200
    }
    # cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2])))
}
ADAM_MLE <- function(X_list, Y_list, Z_list, x_init, alpha = 0.001, beta1 = 0.9, beta2 = 0.999, eta = 10^(-8)) {
  tol <- 10^(-4)
  L <- length(x_init)
  p <- length(x_init) - 2
  gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd)
  }
  m <- rep(0, L)
  v <- rep(0, L)
  x_old <- x_init
  x_new <- x_old - alpha * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list)
  iter <- 1
  while (norm(x_new - x_old, "2") > tol & iter < 5000) {
    x_old <- x_new
    g <- gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list)
    m <- beta1 * m + (1 - beta1) * g
    v <- beta2 * v + (1 - beta2) * (g^2)
    m <- m / (1 - beta1^iter)
    v <- v / (1 - beta2^iter)
    x_new <- x_old - alpha * m / (sqrt(v) + eta)
    iter <- iter + 1
    # cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2])))
}
GD_MLE <- function(X_list, Y_list, Z_list, x_init, step_size) {
  tol <- 10^(-4)
  momentum <- 0.5
  p <- length(x_init) - 2
  gradient <- function(beta, A, B, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + exp(B) * Zi %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- as.numeric(exp(B) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% t(Zi) * exp(B)))
      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd)
  }
  delta_x_old <- 0
  x_old <- x_init
  delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
  x_new <- x_old - delta_x_new
  cnt <- 1
  iter <- 1
  while (norm(x_new - x_old, "2") > tol & iter < 50000) {
    delta_x_old <- delta_x_new
    x_old <- x_new
    delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[p + 2], X_list, Y_list, Z_list) + momentum * delta_x_old
    x_new <- x_old - delta_x_new
    cnt <- cnt + 1
    iter <- iter + 1
    if (cnt > 200) {
      step_size <- 0.8 * step_size
      cnt <- cnt - 200
    }
    #  cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma_sq = exp(x_new[p + 1]), sigma_u_sq = exp(x_new[p + 2])))
}


generate_beta_0 <- function(p = 20, r = 1, sparse_prop = 0.25) {
  rbinom(p, 1, prob = 1 - sparse_prop) * r
}

generate_one_group_X1 <- function(n_i, p, q) {
  design <- t(sapply(1:n_i, function(i) {
    generate_one_X_Z_1(i, p, q)
  }))
  design[, 1:p]
}
generate_one_group_Z1 <- function(n_i, p, q) {
  design <- t(sapply(1:n_i, function(i) {
    generate_one_X_Z_1(i, p, q)
  }))
  design[, (p + 1):(p + q)]
}

generate_one_X_Z_1 <- function(index, p = 10, q = 6) {
  Zij_X <- rnorm(p)
  Zij_Z <- rnorm(q)
  index_p <- (2 * (1:p)) < p
  index_q <- (2 * (1:q)) < q
  X <- Zij_X
  Z <- Zij_Z
  X[index_p] <- 2 * (Zij_X[index_p] > 0) - 1
  Z[index_q] <- 2 * (Zij_Z[index_q] > 0) - 1
  return(c(X, Z))
}
generate_one_X_Z_2 <- function(index, p = 10, q = 6) {
  return(c(rnorm(p), rnorm(q)))
}

generate_one_group_X2 <- function(n_i, p, q) {
  design <- t(sapply(1:n_i, function(i) {
    generate_one_X_Z_2(i, p, q)
  }))
  design[, 1:p]
}
generate_one_group_Z2 <- function(n_i, p, q) {
  design <- t(sapply(1:n_i, function(i) {
    generate_one_X_Z_2(i, p, q)
  }))
  design[, (p + 1):(p + q)]
}
generate_beta_0 <- function(p = 20, r = 1, sparse_prop = 0.25) {
  rbinom(p, 1, prob = 1 - sparse_prop) * r
}

IRS_M_each_people_X <- function(X, M) {
  p <- ncol(X)
  X_syn <- as.matrix(sapply(1:p, function(i) {
    sample(X[, i], M, replace = T)
  }))
  return(X_syn)
}
IRS_M_each_people_Z <- function(Z, M) {
  q <- ncol(Z)
  Z_syn <- as.matrix(sapply(1:q, function(i) {
    sample(Z[, i], M, replace = T)
  }))
  return(Z_syn)
}




GD_ridge_Feb22 <- function(X_list, Y_list, Z_list, x_init, lambda, step_size) {
  tol <- 10^(-5)
  momentum <- 0.5
  p <- ncol(X_list[[1]])
  q <- ncol(Z_list[[1]])
  gradient <- function(beta, A, B_seq, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + Zi %*% diag(exp(B_seq)) %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- sapply(1:length(B_seq), function(index) {
        tep <- rep(0, q)
        tep[index] <- exp(B_seq[index])
        as.numeric(t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% diag(tep) %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% diag(tep) %*% t(Zi)))
      })

      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd + 2 * lambda * c(0, beta[-1], 0, rep(0, q)))
  }
  delta_x_old <- 0
  x_old <- x_init
  delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[(p + 2):length(x_init)], X_list, Y_list, Z_list) + momentum * delta_x_old
  # delta_x_new=step_size*gradient(x_old[1:p],x_old[p+1],x_old[(p+2):length(x_init)],X_list,Y_list,Z_list)+momentum*delta_x_old
  x_new <- x_old - delta_x_new
  cnt <- 1
  iter <- 1
  while (norm(x_new - x_old, "2") > tol & iter < 50000) {
    delta_x_old <- delta_x_new
    x_old <- x_new
    delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[(p + 2):length(x_init)], X_list, Y_list, Z_list) + momentum * delta_x_old
    x_new <- x_old - delta_x_new
    cnt <- cnt + 1
    iter <- iter + 1
    if (cnt > 100) {
      step_size <- 0.8 * step_size
      cnt <- cnt - 100
    }
    # cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma = sqrt(exp(x_new[p + 1])), sigma_u = sqrt(exp(x_new[(p + 2):length(x_init)]))))
}

GD_MLE_Feb22 <- function(X_list, Y_list, Z_list, x_init, step_size) {
  tol <- 10^(-5)
  momentum <- 0.5
  p <- ncol(X_list[[1]])
  q <- ncol(Z_list[[1]])
  gradient <- function(beta, A, B_seq, X_list, Y_list, Z_list) {
    m <- length(X_list)
    n_i <- nrow(X_list[[1]])
    gd <- apply(sapply(1:m, function(i) {
      Xi <- X_list[[i]]
      Yi <- Y_list[[i]]
      Zi <- Z_list[[i]]
      Vi_inv <- ginv(exp(A) * diag(n_i) + Zi %*% diag(exp(B_seq)) %*% t(Zi))

      component1 <- as.numeric(t(Xi) %*% Vi_inv %*% Yi - t(Xi) %*% Vi_inv %*% Xi %*% beta)
      component2 <- as.numeric(exp(A) * t(Yi - Xi %*% beta) %*% Vi_inv %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv * exp(A)))
      component3 <- sapply(1:length(B_seq), function(index) {
        tep <- rep(0, q)
        tep[index] <- exp(B_seq[index])
        as.numeric(t(Yi - Xi %*% beta) %*% Vi_inv %*% Zi %*% diag(tep) %*% t(Zi) %*% Vi_inv %*% (Yi - Xi %*% beta) - tr(Vi_inv %*% Zi %*% diag(tep) %*% t(Zi)))
      })

      return(c(component1, component2, component3))
    }), 1, sum)
    return(-gd)
  }
  delta_x_old <- 0
  x_old <- x_init
  delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[(p + 2):length(x_init)], X_list, Y_list, Z_list) + momentum * delta_x_old
  x_new <- x_old - delta_x_new
  cnt <- 1
  iter <- 1
  while (norm(x_new - x_old, "2") / norm(x_old, "2") > tol & iter < 50000) {
    delta_x_old <- delta_x_new
    x_old <- x_new
    delta_x_new <- step_size * gradient(x_old[1:p], x_old[p + 1], x_old[(p + 2):length(x_init)], X_list, Y_list, Z_list) + momentum * delta_x_old
    x_new <- x_old - delta_x_new
    cnt <- cnt + 1
    iter <- iter + 1
    if (cnt > 100) {
      step_size <- 0.8 * step_size
      cnt <- cnt - 100
    }
    # cat("iteration: ",iter,"\n")
  }
  return(list(beta = x_new[1:p], sigma = sqrt(exp(x_new[p + 1])), sigma_u = sqrt(exp(x_new[(p + 2):length(x_init)]))))
}


function_convenience <- function(m, n_i, p, q, standardize = FALSE, initisGD = FALSE) {
  sigma_u_square <- 1
  sigma_square <- 1
  design_type <- 1
  set.seed(123)
  beta_0 <- generate_beta_0(p, 2, 0.7)
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
  # cat(X_list[[1]][1,1:5],"\n")
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
  V_list <- lapply(1:m, function(index) {
    sigma_square^2 * diag(n_i) + Z_list[[index]] %*% Var_G %*% t(Z_list[[index]])
  })
  random_effect <- mvrnorm(m, rep(0, q), Var_G)
  Y_list <- lapply(1:m, function(index) {
    X_list[[index]] %*% beta_0 + Z_list[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i), diag(n_i))
  })
  mydata <- matrix(0, nc = p + q + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- X_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- Z_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("Y", paste0("X", 1:p), paste0("Z", 1:q), "group")
  GD_MLE_result <- GD_MLE(X_list, Y_list, Z_list, rep(0, p + 1 + 1), 0.01)
  if (initisGD) {
    lmm.MLE <- lmmlasso.default(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
      mydata$group,
      coefInit = list(as.numeric(GD_MLE_result$beta), sqrt(GD_MLE_result$sigma_u), sqrt(GD_MLE_result$sigma_sq)), lambda = 0, startValue = 2,
      pdMat = "pdIdent", CovOpt = "nlm", standardize = F
    )
  } else {
    lmm.MLE <- lmmlasso.default(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
      mydata$group,
      lambda = 0, startValue = 2,
      pdMat = "pdIdent", CovOpt = "nlm", standardize = F
    )
  }

  data.frame(
    MSE_lmm = norm(beta_0 - as.numeric(lmm.MLE$coefficients), "2")^2 / p,
    MSE_GD = norm(beta_0 - as.numeric(GD_MLE_result$beta), "2")^2 / p,
    sigma_u_lmm = lmm.MLE$pars,
    sigma_u_GD = sqrt(GD_MLE_result$sigma_u), sigma_lmm = lmm.MLE$sigma, sigma_GD = sqrt(GD_MLE_result$sigma_sq)
  )
}
function_convenience_ridge <- function(m, n_i, p, q, standardize = FALSE, initisGD = FALSE, lambdainput = 0.1) {
  sigma_u_square <- 1
  sigma_square <- 1
  design_type <- 1
  set.seed(123)
  beta_0 <- generate_beta_0(p, 2, 0.7)
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
  # cat(X_list[[1]][1,1:5],"\n")
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
  V_list <- lapply(1:m, function(index) {
    sigma_square^2 * diag(n_i) + Z_list[[index]] %*% Var_G %*% t(Z_list[[index]])
  })
  random_effect <- mvrnorm(m, rep(0, q), Var_G)
  Y_list <- lapply(1:m, function(index) {
    X_list[[index]] %*% beta_0 + Z_list[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i), diag(n_i))
  })
  mydata <- matrix(0, nc = p + q + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- X_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- Z_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("Y", paste0("X", 1:p), paste0("Z", 1:q), "group")
  GD_MLE_result <- GD_MLE(X_list, Y_list, Z_list, rep(0, p + 1 + 1), 0.01)
  if (initisGD) {
    lmm.MLE <- lmmridge(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
      mydata$group,
      coefInit = list(as.numeric(GD_MLE_result$beta), sqrt(GD_MLE_result$sigma_u), sqrt(GD_MLE_result$sigma_sq)),
      lambda = lambdainput, startValue = 2, pdMat = "pdIdent", CovOpt = "nlm", standardize = F
    )
  } else {
    lmm.MLE <- lmmridge(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
      mydata$group,
      lambda = lambdainput, startValue = 2,
      pdMat = "pdIdent", CovOpt = "nlm", standardize = F
    )
  }

  data.frame(
    MSE_lmm = norm(beta_0 - as.numeric(lmm.MLE$coefficients), "2")^2 / p,
    sigma_u_lmm = lmm.MLE$pars, sigma_lmm = lmm.MLE$sigma
  )
}


















lmm_synthetic_data_generation_submodel <- function(X_list, Y_list, Z_list, sub_model_index_X, sub_model_index_Z, M = NULL, standardize = F) {
  m <- length(X_list)
  n_i <- nrow(X_list[[1]])
  p <- ncol(X_list[[1]])
  q <- ncol(Z_list[[1]])
  X_list_sub <- lapply(1:m, function(i) {
    X_list[[i]][, sub_model_index_X]
  })
  Z_list_sub <- lapply(1:m, function(i) {
    Z_list[[i]][, sub_model_index_Z]
  })
  p_sub <- length(sub_model_index_X)
  q_sub <- length(sub_model_index_Z)
  mydata_sub <- matrix(0, nc = p_sub + q_sub + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata_sub[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata_sub[((i - 1) * n_i + 1):(i * n_i), 2:(p_sub + 1)] <- X_list_sub[[i]]
    mydata_sub[((i - 1) * n_i + 1):(i * n_i), (p_sub + 2):(p_sub + 1 + q_sub)] <- Z_list_sub[[i]]
    mydata_sub[((i - 1) * n_i + 1):(i * n_i), (p_sub + q_sub + 2)] <- i
  }

  lmm.MLE_sub <- lmmlasso.default(as.matrix(mydata_sub[, 2:(p_sub + 1)]), mydata_sub[, 1], as.matrix(mydata_sub[, (p_sub + 2):(p_sub + q_sub + 1)]),
    mydata_sub[, p_sub + q_sub + 2],
    lambda = 0, startValue = 2, pdMat = "pdIdent", CovOpt = "nlm", standardize = standardize
  )
  # print("here")
  beta_0_hat <- as.numeric(lmm.MLE_sub$coefficients)
  Var_G_sub <- lmm.MLE_sub$pars^2 * diag(q_sub)
  sigma_sub <- lmm.MLE_sub$sigma

  #  next we can generate the synthetic X and Z
  if (is.null(M)) M <- round(4 * (p + q) / m)

  X_list_syn <- lapply(1:m, function(index) {
    IRS_M_each_people_X(X_list[[index]], M)
  })
  Z_list_syn <- lapply(1:m, function(index) {
    IRS_M_each_people_Z(Z_list[[index]], M)
  })
  # we can generate synthetic Y

  random_effect_syn <- mvrnorm(m, rep(0, q_sub), Var_G_sub)
  Y_list_syn <- lapply(1:m, function(index) {
    X_list_syn[[index]][, sub_model_index_X] %*% beta_0_hat + Z_list_syn[[index]][, sub_model_index_Z] %*% random_effect_syn[index, ]
    +mvrnorm(1, rep(0, M), sigma_sub^2 * diag(M))
  })

  return(list(X_list_syn = X_list_syn, Y_list_syn = Y_list_syn, Z_list_syn = Z_list_syn, beta_0_hat = beta_0_hat))
}

prediction_given_cv_fitted_object <- function(cv_ridge_result, cv_syn_result, lmm.MLE, X_list_new, Z_list_new) {
  q <- ncol(Z_list_new[[1]])
  random_effect_est_cv_ridge_abs <-
    matrix(cv_ridge_result$random_abs, nc = q, byrow = T)
  random_effect_est_cv_ridge_L2 <-
    matrix(cv_ridge_result$random_L2, nc = q, byrow = T)
  random_effect_est_cv_syn_abs <-
    matrix(cv_syn_result$random_abs, nc = q, byrow = T)
  random_effect_est_cv_syn_L2 <-
    matrix(cv_syn_result$random_L2, nc = q, byrow = T)
  random_effect_est_mle <- matrix(lmm.MLE$random, nc = q, byrow = T)
  betahat_syn_abs <- as.numeric(cv_syn_result$syn_abs$betahat)
  betahat_syn_L2 <- as.numeric(cv_syn_result$syn_L2$betahat)
  betahat_ridge_abs <- as.numeric(cv_ridge_result$ridge_abs$betahat)
  betahat_ridge_L2 <- as.numeric(cv_ridge_result$ridge_L2$betahat)
  betahat_MLE <- as.numeric(lmm.MLE$coefficients)
  N <- length(X_list_new)
  predGrp_MLE <- list()
  length(predGrp_MLE) <- N
  predGrp_syn_abs <- list()
  length(predGrp_syn_abs) <- N
  predGrp_syn_L2 <- list()
  length(predGrp_syn_L2) <- N
  predGrp_ridge_abs <- list()
  length(predGrp_ridge_abs) <- N
  predGrp_ridge_L2 <- list()
  length(predGrp_ridge_L2) <- N
  for (i in 1:N)
  {
    predGrp_MLE[[i]] <-
      X_list_new[[i]] %*% betahat_MLE + Z_list_new[[i]] %*% random_effect_est_mle[i, ]
    predGrp_syn_abs[[i]] <-
      X_list_new[[i]] %*% betahat_syn_abs + Z_list_new[[i]] %*% random_effect_est_cv_syn_abs[i, ]
    predGrp_syn_L2[[i]] <-
      X_list_new[[i]] %*% betahat_syn_L2 + Z_list_new[[i]] %*% random_effect_est_cv_syn_L2[i, ]
    predGrp_ridge_abs[[i]] <-
      X_list_new[[i]] %*% betahat_ridge_abs + Z_list_new[[i]] %*% random_effect_est_cv_ridge_abs[i, ]
    predGrp_ridge_L2[[i]] <-
      X_list_new[[i]] %*% betahat_ridge_L2 + Z_list_new[[i]] %*% random_effect_est_cv_ridge_L2[i, ]
  }
  return(list(
    predGrp_MLE = predGrp_MLE, predGrp_syn_abs = predGrp_syn_abs, predGrp_syn_L2 = predGrp_syn_L2,
    predGrp_ridge_abs = predGrp_ridge_abs, predGrp_ridge_L2 = predGrp_ridge_L2
  ))
}


# this function is used to evaluate the prediction accuracy for MLE, pred_syn_abs,pred_syn_syn_L2,pred_ridge_abs, pred_ridge_L2
pred_error_based_on_Y_new <- function(predGrp_MLE, predGrp_syn_abs, predGrp_syn_L2, predGrp_ridge_abs, predGrp_ridge_L2, Y_new_list) {
  N <- length(Y_new_list)
  residGrp_MLE <- list()
  length(residGrp_MLE) <- N
  residGrp_syn_abs <- list()
  length(residGrp_syn_abs) <- N
  residGrp_syn_L2 <- list()
  length(residGrp_syn_L2) <- N
  residGrp_ridge_abs <- list()
  length(residGrp_ridge_abs) <- N
  residGrp_ridge_L2 <- list()
  length(residGrp_ridge_L2) <- N
  for (i in 1:N) {
    residGrp_MLE[[i]] <- predGrp_MLE[[i]] - Y_new_list[[i]]

    residGrp_syn_abs[[i]] <- predGrp_syn_abs[[i]] - Y_new_list[[i]]

    residGrp_syn_L2[[i]] <- predGrp_syn_L2[[i]] - Y_new_list[[i]]

    residGrp_ridge_abs[[i]] <- predGrp_ridge_abs[[i]] - Y_new_list[[i]]

    residGrp_ridge_L2[[i]] <- predGrp_ridge_L2[[i]] - Y_new_list[[i]]
  }
  return(c(
    pred_error_MLE_abs = sum(abs(unlist(residGrp_MLE))),
    pred_error_MLE_L2 = sum(unlist(residGrp_MLE)^2),
    pred_error_syn_abs = sum(abs(unlist(residGrp_syn_abs))),
    pred_error_syn_L2 = sum(unlist(residGrp_syn_L2)^2),
    pred_error_ridge_abs = sum(abs(unlist(residGrp_ridge_abs))),
    pred_error_ridge_L2 = sum(unlist(residGrp_ridge_L2)^2)
  ))
}
MSE_given_cv_fitted_object <- function(cv_ridge_result, cv_syn_result, lmm.MLE, beta_0) {
  betahat_syn_abs <- as.numeric(cv_syn_result$syn_abs$betahat)
  betahat_syn_L2 <- as.numeric(cv_syn_result$syn_L2$betahat)
  betahat_ridge_abs <- as.numeric(cv_ridge_result$ridge_abs$betahat)
  betahat_ridge_L2 <- as.numeric(cv_ridge_result$ridge_L2$betahat)
  betahat_MLE <- as.numeric(lmm.MLE$coefficients)
  p <- length(beta_0)
  return(c(
    MSE_MLE_abs = sum(abs(beta_0 - betahat_MLE)) / p,
    MSE_MLE_L2 = sum((beta_0 - betahat_MLE)^2) / p,
    MSE_syn_abs = sum(abs(beta_0 - betahat_syn_abs)) / p,
    MSE_syn_L2 = sum((beta_0 - betahat_syn_L2)^2) / p,
    MSE_ridge_abs = sum(abs(beta_0 - betahat_ridge_abs)) / p,
    MSE_ridge_L2 = sum((beta_0 - betahat_ridge_L2)^2) / p
  ))
}


Apr12_one_time_simulation <- function(m, n_i, p, q, sigma_u_square, sigma_square, design_type = 1, standardize = FALSE) {
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
  V_list <- lapply(1:m, function(index) {
    sigma_square^2 * diag(n_i) + Z_list[[index]] %*% Var_G %*% t(Z_list[[index]])
  })
  random_effect <- mvrnorm(m, rep(0, q), Var_G)
  Y_list <- lapply(1:m, function(index) {
    X_list[[index]] %*% beta_0 + Z_list[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i), diag(n_i))
  })

  mydata <- matrix(0, nc = p + q + 2, nr = m * n_i)
  for (i in 1:m) {
    mydata[((i - 1) * n_i + 1):(i * n_i), 1] <- Y_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- X_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- Z_list[[i]]
    mydata[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("Y", paste0("X", 1:p), paste0("Z", 1:q), "group")
  lmm.MLE <- lmmlasso.default(as.matrix(mydata[, 2:(p + 1)]), mydata$Y, as.matrix(mydata[, (p + 2):(p + q + 1)]),
    mydata$group,
    lambda = 0, startValue = 2, pdMat = "pdIdent", CovOpt = "nlm", standardize = F
  )

  my_syndata <- lmm_synthetic_data_generation_submodel(X_list, Y_list, Z_list, 1:2, 1:2)
  tau_seq <- 10^seq(-6, 0, by = 0.5)
  lambda_seq <- tau_seq
  whole_set <- five_fold_CV_separation(X_list, Z_list, Y_list)
  cv_ridge_result <- five_fold_CV_lmmridge(X_list, Z_list, Y_list, whole_set, lambda_seq)
  cv_syn_result <- five_fold_CV_lmmsyn(X_list, Z_list, Y_list, whole_set, tau_seq, my_syndata)
  n_i_new <- 5
  X_list_new <- lapply(1:m, function(index) {
    if (design_type == 1) {
      x <- cbind(1, generate_one_group_X1(n_i_new, p - 1, q))
      if (standardize) {
        meanx <- apply(x[, -1], 2, mean)
        sdx <- apply(x[, -1], 2, sd)
        x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
      }
      return(x)
    } else {
      x <- cbind(1, generate_one_group_X2(n_i_new, p - 1, q))
      if (standardize) {
        meanx <- apply(x[, -1], 2, mean)
        sdx <- apply(x[, -1], 2, sd)
        x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
      }
      return(x)
    }
  })
  Z_list_new <- lapply(1:m, function(index) {
    if (design_type == 1) {
      x <- generate_one_group_Z1(n_i_new, p, q)
      if (standardize) {
        meanx <- apply(x[, -1], 2, mean)
        sdx <- apply(x[, -1], 2, sd)
        x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
      }
      return(x)
    } else {
      x <- generate_one_group_Z2(n_i_new, p, q)
      if (standardize) {
        meanx <- apply(x[, -1], 2, mean)
        sdx <- apply(x[, -1], 2, sd)
        x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
      }
      return(x)
    }
  })
  Y_list_new <- lapply(1:m, function(index) {
    X_list_new[[index]] %*% beta_0 + Z_list_new[[index]] %*% random_effect[index, ] + mvrnorm(1, rep(0, n_i_new), diag(n_i_new))
  })

  prediction_result <- prediction_given_cv_fitted_object(cv_ridge_result, cv_syn_result, lmm.MLE, X_list_new, Z_list_new)
  pred_error_table <- pred_error_based_on_Y_new(
    prediction_result$predGrp_MLE, prediction_result$predGrp_syn_abs, prediction_result$predGrp_syn_L2,
    prediction_result$predGrp_ridge_abs, prediction_result$predGrp_ridge_L2, Y_list_new
  )
  MSE_table <- MSE_given_cv_fitted_object(cv_ridge_result, cv_syn_result, lmm.MLE, beta_0)
  rere <- list(pred_error_table, MSE_table)
  save(rere, file = paste0("error_table_P", p))

  return(rere)
}






MAP_syn_iteration_with_sigma_alpha <- function(X_list, Y_list, Z_list, X_list_syn, Y_list_syn, Z_list_syn, prior_alpha0 = 1, prior_beta0 = 1,
                                               prior_alpha1 = 1, prior_beta1 = 1, tau, beta_init, sigma_init = NULL, sigma_alpha_init = NULL, grid_sigma, MAX_ITER = 1000) {
  tol <- 1e-06
  p <- ncol(X_list[[1]])
  n_i <- nrow(X_list[[1]])
  M <- nrow(X_list_syn[[1]])
  q <- ncol(Z_list[[1]])
  m <- length(X_list)

  g_former_part <- rep(0, p)
  G_former_part <- matrix(0, p, p)

  for (i in 1:m) {
    g_former_part <- g_former_part + t(X_list[[i]]) %*% Y_list[[i]] + tau / M * t(X_list_syn[[i]]) %*% Y_list_syn[[i]]
    G_former_part <- G_former_part + t(X_list[[i]]) %*% X_list[[i]] + tau / M * t(X_list_syn[[i]]) %*% X_list_syn[[i]]
  }
  c_list <- lapply(1:m, function(index) {
    -t(Z_list[[index]]) %*% X_list[[index]] - tau / M * t(Z_list_syn[[index]]) %*% X_list_syn[[index]]
  })
  d_list <- lapply(1:m, function(index) {
    t(Z_list[[index]]) %*% Y_list[[index]] + tau / M * t(Z_list_syn[[index]]) %*% Y_list_syn[[index]]
  })
  ZTZ_list <- lapply(1:m, function(index) {
    t(Z_list[[index]]) %*% Z_list[[index]]
  })
  ZTZ.star_list <- lapply(1:m, function(index) {
    t(Z_list_syn[[index]]) %*% Z_list_syn[[index]]
  })

  constant_sigma_power <- -n_i * m - n_i * tau * m / M + q * m - 2 * (prior_alpha0 + 1)
  constant_sigma_alpha_power <- -q * m - 2 * (prior_alpha1 + 1)

  each_sigma_log_posterior2 <- function(sigma_cur, beta_cur, sigma_alpha_cur) {
    log_sigma_cur <- log(sigma_cur)
    bi_list_beta_cur <- lapply(1:m, function(index) {
      t(Z_list[[index]]) %*% (Y_list[[index]] - X_list[[index]] %*% beta_cur) + tau / M * t(Z_list_syn[[index]]) %*% (Y_list_syn[[index]] - X_list_syn[[index]] %*% beta_cur)
    })

    RSS_long_term_cur_with_det <- sum(sapply(1:m, function(index) {
      bi <- bi_list_beta_cur[[index]]
      Bi <- ZTZ_list[[index]] + tau / M * ZTZ.star_list[[index]] + sigma_cur^2 / sigma_alpha_cur^2 * diag(q)
      tem <- prior_beta0 + norm(Y_list[[index]] - X_list[[index]] %*% beta_cur, "2")^2 + tau / M * norm(Y_list_syn[[index]] - X_list_syn[[index]] %*% beta_cur, "2")^2 - quad.form.inv(Bi, bi)
      -0.5 * tem / sigma_cur^2 + 0.5 * as.numeric(determinant(ginv(Bi))$mod)
    }))
    constant_sigma_power * log_sigma_cur + RSS_long_term_cur_with_det
  }

  each_sigma_alpha_log_posterior2 <- function(sigma_cur, beta_cur, sigma_alpha_cur) {
    log_sigma_alpha_cur <- log(sigma_alpha_cur)
    bi_list_beta_cur <- lapply(1:m, function(index) {
      t(Z_list[[index]]) %*% (Y_list[[index]] - X_list[[index]] %*% beta_cur) + tau / M * t(Z_list_syn[[index]]) %*% (Y_list_syn[[index]] - X_list_syn[[index]] %*% beta_cur)
    })

    RSS_long_term_cur_with_det <- sum(sapply(1:m, function(index) {
      bi <- bi_list_beta_cur[[index]]
      Bi <- ZTZ_list[[index]] + tau / M * ZTZ.star_list[[index]] + sigma_cur^2 / sigma_alpha_cur^2 * diag(q)
      tem <- -quad.form.inv(Bi, bi)
      -0.5 * tem / sigma_cur^2 + 0.5 * as.numeric(determinant(ginv(Bi))$mod)
    }))
    constant_sigma_alpha_power * log_sigma_alpha_cur + RSS_long_term_cur_with_det - prior_beta1 / 2 / sigma_alpha_cur^2
  }


  beta_new <- beta_init
  sigma_new <- sigma_init

  grid_optimize_result_sigma_alpha <- sapply(grid_sigma, function(sig) {
    each_sigma_alpha_log_posterior2(sigma_new, beta_new, sig)
  })
  sigma_alpha_new <- grid_sigma[which.max(grid_optimize_result_sigma_alpha)]

  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    iter <- iter + 1
    beta_old <- beta_new
    sigma_old <- sigma_new
    sigma_alpha_old <- sigma_alpha_new

    g <- g_former_part
    G <- G_former_part
    for (i in 1:m) {
      Bi <- ZTZ_list[[i]] + tau / M * ZTZ.star_list[[i]] + sigma_old^2 / sigma_alpha_old^2 * diag(q)
      g <- g + quad.3form.inv(Bi, c_list[[i]], d_list[[i]])
      G <- G - quad.form.inv(Bi, c_list[[i]])
    }
    G_inv <- ginv(G)
    beta_new <- as.numeric(G_inv %*% g)

    grid_optimize_result_sigma <- sapply(grid_sigma, function(sig) {
      each_sigma_log_posterior2(sig, beta_new, sigma_alpha_old)
    })
    sigma_new <- grid_sigma[which.max(grid_optimize_result_sigma)]

    grid_optimize_result_sigma_alpha <- sapply(grid_sigma, function(sig) {
      each_sigma_alpha_log_posterior2(sigma_new, beta_new, sig)
    })
    sigma_alpha_new <- grid_sigma[which.max(grid_optimize_result_sigma_alpha)]




    if (norm(c(beta_old, sigma_old, sigma_alpha_old) - c(beta_new, sigma_new, sigma_alpha_new), "2") > tol * norm(c(beta_old, sigma_old, sigma_alpha_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  # print(iter)
  if (iter > (MAX_ITER - 2)) {
    print("not converge")
  }
  return(list(betahat = beta_new, sigma_hat = sigma_new, sigma_alpha_hat = sigma_alpha_new))
}

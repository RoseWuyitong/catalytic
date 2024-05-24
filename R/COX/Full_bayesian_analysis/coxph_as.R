# copy the theorem from the website: Newton-Raphson https://medium.com/analytics-vidhya/implementing-the-cox-model-in-r-b1292d6ab6d2
library(survival)
# library(ggplot2)
library(dplyr)
library(glmnet)
library(emulator)
ITERS <- 25
INITIAL_BETA <- 0
# Get the risk function
GetRiskSet <- function(time_of_interest,
                       entry_vector,
                       time_vector,
                       event_vector) {
  return(which((time_of_interest >= entry_vector) & ((time_vector == time_of_interest & event_vector == 1) | (time_vector + 1e-08 > time_of_interest))))
}
# compute the gradient of the likelihood
CoxGradient <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  p <- ncol(Xs)

  gradient <- apply(
    cbind(Ts, Xs), 1,
    function(df) {
      df <- matrix(df, nrow = 1)
      ts <- df[, 1]
      xs <- df[, 2:ncol(df)]
      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event), ] %>%
        matrix(ncol = ncol(Xs))

      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
      t2 <- sum(exp(X_risk_set %*% beta))

      return(xs - t1 / t2)
    }
  ) %*% event # %>%
  # matrix(nrow = p) %>%
  # rowSums()
  # print(dim(gradient));stop()
  return(gradient)
}

# compute the variance
CoxVariance <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  p <- ncol(Xs)

  variance <- apply(
    cbind(Ts, Xs), 1,
    function(df) {
      df <- matrix(df, nrow = 1)
      ts <- df[, 1]
      xs <- df[, 2:ncol(df)]
      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event), ] %>%
        matrix(ncol = ncol(Xs))

      sum1 <- sum(exp(X_risk_set %*% beta))
      sum2 <- rowSums(t(X_risk_set^2) %*% exp(X_risk_set %*% beta))

      sum3 <- rowSums(t(X_risk_set) %*% exp(X_risk_set %*% beta))^2
      # print(sum2);stop()
      return(-(sum1 * sum2 - sum3) / sum1^2)
    }
  ) %>%
    matrix(nrow = p) %>%
    rowSums()
  return(-1 / variance)
}
# plot the gradient
GradientPlot <- function(gradients) {
  gradient_track <- reshape2::melt(gradients) %>%
    `names<-`(c("iteration", "beta", "gradient"))

  p <- ggplot(gradient_track, aes(y = gradient, x = iteration)) +
    geom_line() +
    facet_wrap(~ factor(beta))

  return(p)
}

# plot and make comparison with standard package
BetaPlot <- function(store_betas, store_variance, coxph_output) {
  variance_track <- reshape2::melt(store_variance) %>%
    `names<-`(c("iteration", "variable", "coxph_AS_variance"))

  plot_df <- reshape2::melt(store_betas) %>%
    `names<-`(c("iteration", "variable", "coxph_AS_beta")) %>%
    left_join(coxph_output, by = "variable") %>%
    left_join(variance_track, by = c("variable", "iteration")) %>%
    mutate(
      coxph_AS_lci = coxph_AS_beta - 1.96 * sqrt(coxph_AS_variance),
      coxph_AS_uci = coxph_AS_beta + 1.96 * sqrt(coxph_AS_variance)
    )

  ggplot(plot_df) +
    geom_line(aes(y = coxph_AS_beta, x = iteration), size = .8) +
    geom_line(aes(y = coxph_AS_lci, x = iteration), linetype = "dashed", size = .8) +
    geom_line(aes(y = coxph_AS_uci, x = iteration), linetype = "dashed", size = .8) +
    geom_line(aes(y = coxph_beta, x = iteration, color = "red")) +
    geom_line(aes(y = coxph_lci, x = iteration), linetype = "dashed", color = "red") +
    geom_line(aes(y = coxph_uci, x = iteration), linetype = "dashed", color = "red") +
    facet_wrap(~ factor(variable_name)) +
    ylab("beta estimate") +
    theme(legend.position = "none")
}
# our own implementation of coxph, combine function above
coxph_AS <- function(formula, data) {
  # Fit using coxph()

  model <- coxph(formula, data)

  cox_beta <- data.frame(
    variable = 1:length(coef(model)),
    coxph_beta = coef(model)
  )
  cox_ci <- cbind(1:length(coef(model)), confint(model)) %>%
    `colnames<-`(c("variable", "coxph_lci", "coxph_uci")) %>%
    data.frame()
  coxph_output <- left_join(cox_beta, cox_ci, by = "variable") %>%
    mutate(variable_name = names(coef(model)))

  # Load data
  x <- data %>%
    select(all.vars(formula[-1])) %>%
    as.matrix()
  lhs <- setdiff(all.vars(formula), all.vars(formula[-1]))

  if (length(lhs) == 1) { # no censoring
    times <- data %>%
      select(all.vars(formula)[1]) %>%
      as.matrix()
    event <- rep(1, dim(times)[1]) %>% as.matrix()
    entry <- rep(0, dim(times)[1]) %>% as.matrix()
  } else if (length(lhs) == 2) { # censoring but no entry time specified
    times <- data %>%
      select(all.vars(formula)[1]) %>%
      as.matrix()
    event <- data %>%
      select(all.vars(formula)[2]) %>%
      as.matrix()
    entry <- rep(0, dim(times)[1]) %>% as.matrix()
  } else if (length(lhs) == 3) { # censoring but no entry time specified
    times <- data %>%
      select(all.vars(formula)[2]) %>%
      as.matrix()
    event <- data %>%
      select(all.vars(formula)[3]) %>%
      as.matrix()
    entry <- data %>%
      select(all.vars(formula)[1]) %>%
      as.matrix()
  }
  # Initialize matrices
  store_gradient <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  store_betas <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  store_variance <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  beta <- matrix(rep(INITIAL_BETA, ncol(x)), nrow = ncol(x))

  # Newton-Raphson iterations

  for (i in 1:ITERS) {
    store_gradient[i, ] <- CoxGradient(beta, x, entry, times, event)
    store_variance[i, ] <- CoxVariance(beta, x, entry, times, event)
    store_betas[i, ] <- beta <- beta + store_gradient[i, ] * store_variance[i, ]
  }

  # Plot
  beta_plot <- BetaPlot(store_betas, store_variance, coxph_output)
  gradient_plot <- GradientPlot(store_gradient)

  # Final output data.frame
  coxph_output$coxph_AS_beta <- store_betas[i, ]
  coxph_output$coxph_AS_lci <- store_betas[i, ] - 1.96 * sqrt(store_variance[i, ])
  coxph_output$coxph_AS_uci <- store_betas[i, ] + 1.96 * sqrt(store_variance[i, ])

  coxph_output <- coxph_output %>%
    transmute(
      variable = variable_name,
      beta = coxph_beta,
      beta_AS = coxph_AS_beta,
      LCI = coxph_lci,
      LCI_AS = coxph_AS_lci,
      UCI = coxph_uci,
      UCI_AS = coxph_AS_uci
    )

  return(list(
    beta_plot = beta_plot,
    gradient_plot = gradient_plot,
    model_output = coxph_output
  ))
}


CoxGradient <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  p <- ncol(Xs)

  gradient <- apply(
    cbind(Ts, Xs), 1,
    function(df) {
      df <- matrix(df, nrow = 1)
      ts <- df[, 1]
      xs <- df[, 2:ncol(df)]
      risk_index <- GetRiskSet(ts, entry, Ts, event)
      # print(risk_index)
      if (length(risk_index) > 1) {
        X_risk_set <- Xs[risk_index, ]
      } else {
        X_risk_set <- matrix(Xs[risk_index, ], nc = ncol(Xs))
      }
      # X_risk_set <-   ifelse(length(risk_index)>1,Xs[risk_index,],matrix(Xs[risk_index,],nc = ncol(Xs))) #Xs[risk_index,] %>% #
      # matrix(ncol = ncol(Xs))
      # print(length(risk_index)>1);print(Xs[risk_index,]);print(X_risk_set);stop()
      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
      # print(t1);stop()
      t2 <- sum(exp(X_risk_set %*% beta))

      return(xs - t1 / t2)
    }
  ) %*% event
  return(gradient) # this gradient can also be used into GD
}
CoxHessian <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  hessian <- matrix(0, p, p)
  uncersor_set <- which(as.numeric(event) == 1)
  for (i in uncersor_set) {
    ts <- Ts[i]
    xs <- Xs[i, ]
    riskindex <- GetRiskSet(ts, entry, Ts, event)
    if (length(riskindex) > 1) {
      X_risk_set <- Xs[riskindex, ]
    } else {
      X_risk_set <- matrix(Xs[riskindex, ], nc = ncol(Xs))
    }
    # X_risk_set <- Xs[riskindex,] %>%
    #  matrix(ncol = ncol(Xs))
    theta <- exp(X_risk_set %*% beta)
    t1 <- t(X_risk_set) %*% theta
    t2 <- sum(theta)

    xxxt <- t(X_risk_set)
    eee <- sapply(1:length(riskindex), function(i) (xxxt[, i] * sqrt(theta[i])))
    hessian <- hessian + crossprod(t(eee)) / t2 - t1 %*% t(t1) / t2^2
  }

  return(-hessian) # slow is acceptable
}
CoxPL <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  entry <- rep(0, length(Ts))
  uncersor_set <- which(as.numeric(event) == 1)
  each_PL <- sapply(uncersor_set, function(i) {
    ts <- Ts[i]
    xs <- Xs[i, ]
    riskindex <- GetRiskSet(ts, entry, Ts, event)
    if (length(riskindex) > 1) {
      X_risk_set <- Xs[riskindex, ]
    } else {
      X_risk_set <- matrix(Xs[riskindex, ], nc = ncol(Xs))
    }
    # X_risk_set <- Xs[riskindex,] %>%
    #  matrix(ncol = ncol(Xs))
    theta <- exp(X_risk_set %*% beta)
    return(sum(xs * beta) - log(sum(theta)))
  })
  return(sum(each_PL))
}




newton_cox <- function(tol = 1e-7, MAX_ITER = 25, initial_values, TT, status, X) {
  n <- length(TT)
  entry <- rep(0, n)
  beta_old <- initial_values
  H_old <- CoxHessian(beta_old, X, entry, TT, status)
  gd_old <- CoxGradient(beta_old, X, entry, TT, status)
  beta_new <- beta_old - qr_solve_my(H_old, gd_old)
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    beta_old <- beta_new
    H_old <- CoxHessian(beta_old, X, entry, TT, status)
    gd_old <- CoxGradient(beta_old, X, entry, TT, status)
    beta_new <- beta_old - qr_solve_my(H_old, gd_old)
    iter <- iter + 1
    if (sum(is.na(beta_new)) >= 1) {
      return(initial_values)
    }
    if (norm(as.numeric(beta_new), "2") > 10000) {
      return(initial_values)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) cat("not converge \n")
  if (iter < MAX_ITER) cat("number of iteration to convergence: ", iter, "\n")
  return(as.numeric(beta_new))
}



mimic_glmnet_cox_MLE_syn <- function(TT, status, X, my_syndata, init_beta, tau, tol = 1e-06) {
  TT.star <- my_syndata$T.syn
  X.star <- my_syndata$X.syn
  status.star <- my_syndata$status.syn
  M <- length(TT.star)

  p <- ncol(X)
  beta_old <- rnorm(p) + 1
  beta_new <- init_beta
  # beta_new = rep(0,p)
  y <- Surv(TT, status)
  y.star <- Surv(TT.star, status.star)
  while (norm(beta_new - beta_old, "2") > tol * norm(beta_old, "2")) {
    # print("enter?")
    beta_old <- beta_new
    eta <- X %*% beta_old
    eta.star <- X.star %*% beta_old
    derivative <- coxgrad(eta, y, diag.hessian = T)
    derivative.star <- coxgrad(eta.star, y.star, diag.hessian = T)
    W <- attr(derivative, "diag")
    W.star <- attr(derivative.star, "diag")
    W[which(W == 0)] <- mean(W[which(W != 0)])
    W.star[which(W.star == 0)] <- mean(W.star[which(W.star != 0)])
    Z <- as.numeric(eta - derivative / W)
    Z.star <- as.numeric(eta.star - derivative.star / W.star)
    beta_new <- qr_solve_my(
      quad.form(diag(W), X) + tau / M * quad.form(diag(W.star), X.star),
      quad.3form(diag(W), X, Z) + tau / M * quad.3form(diag(W.star), X.star, Z.star)
    )
    # print(W)
    if (norm(beta_new, "2") > 1000 | sum(is.na(beta_new)) > 0) stop("error in the while loop, explode!")
  }
  as.numeric(beta_new)
}

mimic_glmnet_cox_MLE_syn_apr18 <- function(TT, status, X, my_syndata, tau, tol = 1e-06) {
  TT.star <- my_syndata$T.syn
  X.star <- my_syndata$X.syn
  status.star <- my_syndata$status.syn
  TT.tilde <- c(TT, TT.star)
  X.tilde <- rbind(X, X.star)
  status.tilde <- c(status, status.star)

  p <- ncol(X)
  beta_old <- rnorm(p)
  beta_new <- rep(0, p)
  y.tilde <- Surv(TT.tilde, status.tilde)
  weight_tau <- c(rep(1, nrow(X)), rep(tau, nrow(X.star)))

  while (norm(beta_new - beta_old, "2") > tol * norm(beta_old, "2")) {
    beta_old <- beta_new
    eta.tilde <- X.tilde %*% beta_old
    derivative.tilde <- coxgrad(eta.tilde, y.tilde, diag.hessian = T)

    W <- attr(derivative.tilde, "diag")
    W.tilde <- W * weight_tau
    Z.tilde <- as.numeric(eta.tilde - derivative.tilde / W)
    beta_new <- qr_solve_my(quad.form(diag(W.tilde), X.tilde), quad.3form(diag(W.tilde), X.tilde, Z.tilde))
    if (norm(beta_new, "2") > 1000 | sum(is.na(beta_new)) > 0) stop("error in the while loop, explode!")
  }
  as.numeric(beta_new)
}

mimic_glmnet_cox_MLE <- function(TT, status, X, tol = 1e-06) {
  p <- ncol(X)
  beta_old <- rnorm(p)
  beta_new <- rep(0, p)
  y <- Surv(TT, status)
  while (norm(beta_new - beta_old, "2") > tol * norm(beta_old, "2")) {
    beta_old <- beta_new
    eta <- X %*% beta_old
    derivative <- coxgrad(eta, y, diag.hessian = T)
    W <- attr(derivative, "diag")
    Z <- as.numeric(eta - derivative / W)
    beta_new <- qr_solve_my(quad.form(diag(W), X), quad.3form(diag(W), X, Z))
    if (norm(beta_new, "2") > 1000 | sum(is.na(beta_new)) > 0) stop("error in the while loop, explode!")
  }
  as.numeric(beta_new)
}




IRS_M_each_people_X <- function(X, M) {
  p <- ncol(X)
  X_syn <- as.matrix(sapply(1:p, function(i) {
    sample(X[, i], M, replace = T)
  }))
  return(X_syn)
}
# write a function to generate synthetic data, beta_0_hat can be from sub model
# we will not generate censored data first, and T.syn only depend on uncensored data
generate_synthetic_data <- function(M = 200, sub_model_index = seq(1, 5), TT, X, status, local_cnt = 5) {
  X.syn <- IRS_M_each_people_X(X, M)
  fit <- glmnet(X[, sub_model_index], Surv(TT, status), family = "cox", alpha = 0, lambda = 0, standardize = F)
  beta_hat <- as.numeric(coef(fit))
  d_star <- sapply(1:M, function(i) {
    sum(X.syn[i, ][sub_model_index] * beta_hat)
  })
  d <- sapply(1:nrow(X), function(i) {
    sum(X[i, ][sub_model_index] * beta_hat)
  })
  d_event <- d[which(status == 1)]
  T_event <- TT[which(status == 1)]
  # print(d_event)
  T.syn <- rep(0, M)
  # synthetic response is local average instead of global weighted average since some of value is big
  for (i in 1:M) {
    whole_dist <- sum(sapply(1:sum(status), function(j) {
      abs(d_star[i] - d_event[j])
    }))
    pairdist <- sapply(1:sum(status), function(j) {
      abs(d_star[i] - d_event[j])
    })
    T.syn[i] <- mean(T_event[order(pairdist)[1:local_cnt]])
  }
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}



compare_cox_oracle <- function(TT, status, X, lambda_seq, tau_seq, my_syndata, beta_0) {
  p <- ncol(X)
  if (length(tau_seq) != length(lambda_seq)) stop("lambda_seq should have same length with tau_seq")
  lambda_error <- rep(0, length(tau_seq))
  tau_error <- rep(0, length(tau_seq))
  # tau_seq=lambda_seq for simplicity
  cv_error_each_tau_lambda <- function(tau) {
    fit.ridge_cur_lambda <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = tau, standardize = F)
    beta_syn_this_tau <- newton_cox_syn(tol = 1e-7, MAX_ITER = 25, initial_values = rep(0, ncol(X)), TT, status, X, my_syndata, tau)
    beta_ridge_this_lambda <- as.numeric(coef(fit.ridge_cur_lambda))
    # print(c(norm(beta_0-beta_syn_this_tau,"2")^2/p,norm(beta_0-beta_ridge_this_lambda,"2")^2/p))
    c(norm(beta_0 - beta_syn_this_tau, "2")^2 / p, norm(beta_0 - beta_ridge_this_lambda, "2")^2 / p)
  }
  tau_lambda_error <- sapply(tau_seq, cv_error_each_tau_lambda)
  tau_error <- tau_lambda_error[1, ]
  lambda_error <- tau_lambda_error[2, ]

  tau_min <- tau_seq[which.min(tau_error)]
  lambda_min <- lambda_seq[which.min(lambda_error)]
  cat("best tau and lambda are : ", tau_min, lambda_min, "\n")
  # MLE result
  fit_ML <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = 0, standardize = F)
  beta_ML <- as.numeric(coef(fit_ML))
  # print(dim(X));print(TT);print(status)
  # ridge result
  fit_ridge <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = lambda_min, standardize = F)
  beta_ridge <- as.numeric(coef(fit_ridge))
  # catlaytic result
  beta_syn <- newton_cox_syn(tol = 1e-7, MAX_ITER = 25, initial_values = rep(0, ncol(X)), TT, status, X, my_syndata, tau_min)
  # beta_table=data.frame(beta_truth=beta_0,beta_ML=beta_ML,beta_ridge=beta_ridge,beta_syn=beta_syn)
  error_table <- c(ML.MSE = norm(beta_ML - beta_0, "2")^2 / p, ridge.MSE = norm(beta_ridge - beta_0, "2")^2 / p, syn.MSE = norm(beta_syn - beta_0, "2")^2 / p)
  list(error_table = error_table)
}


compare_cox_cv <- function(TT, status, X, lambda_seq, tau_seq, nfold = 10, my_syndata, beta_0) {
  p <- ncol(X)
  cv.ridge <- cv.glmnet(X, Surv(TT, status),
    family = "cox", alpha = 0,
    lambda = lambda_seq, standardize = F, nfolds = nfold
  )
  # print("here")
  tau_error <- rep(0, length(tau_seq))
  cv_error_each_tau <- function(tau) {
    validation_error <- rep(0, nfold)
    validation_index_list <- split(1:length(TT), rep_len(1:nfold, length(TT)))
    for (i in 1:nfold) {
      test_index <- validation_index_list[[i]]
      beta_syn_this_fold <- newton_cox_syn(
        tol = 1e-7, MAX_ITER = 25,
        initial_values = rep(0, ncol(X)), TT[-test_index], status[-test_index], X[-test_index, ], my_syndata, tau
      )
      validation_error[i] <- coxnet.deviance(y = Surv(TT[test_index], status[test_index]), x = X[test_index, ], beta = as.numeric(beta_syn_this_fold))
    }
    mean(validation_error)
  }
  tau_error <- sapply(tau_seq, cv_error_each_tau)
  tau_min <- tau_seq[which.min(tau_error)]

  # MLE result
  fit_ML <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = 0, standardize = F)
  beta_ML <- as.numeric(fit_ML$coefficients)
  # ridge result
  beta_ridge <- as.numeric(coef(cv.ridge))
  # catlaytic result
  beta_syn <- newton_cox_syn(tol = 1e-7, MAX_ITER = 25, initial_values = rep(0, ncol(X)), TT, status, X, my_syndata, tau_min)
  beta_table <- data.frame(beta_truth = beta_0, beta_ML = beta_ML, beta_ridge = beta_ridge, beta_syn = beta_syn)
  error_table <- c(ML.MSE = norm(beta_ML - beta_0, "2")^2 / p, ridge.MSE = norm(beta_ridge - beta_0, "2")^2 / p, syn.MSE = norm(beta_syn - beta_0, "2")^2 / p)
  list(beta_table = beta_table, error_table = error_table)
}

compare_glmnet_cox_oracle_mar29 <- function(TT, status, X, lambda_seq, tau_seq, my_syndata, beta_0) {
  p <- ncol(X)
  lambda_error <- rep(0, length(lambda_seq))
  tau_error <- rep(0, length(tau_seq))
  # tau_seq=lambda_seq for simplicity
  cv_error_each_tau <- function(tau) {
    beta_syn_this_tau <- mimic_glmnet_cox_MLE_syn(TT, status, X, my_syndata, tau)
    norm(beta_0 - beta_syn_this_tau, "2")^2 / p
  }
  cv_error_each_lambda <- function(lam) {
    fit.ridge_cur_lambda <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = lam, standardize = F)
    beta_ridge_this_lambda <- as.numeric(coef(fit.ridge_cur_lambda))
    norm(beta_0 - beta_ridge_this_lambda, "2")^2 / p
  }
  tau_error <- sapply(tau_seq, cv_error_each_tau)
  lambda_error <- sapply(lambda_seq, cv_error_each_lambda)


  tau_min <- tau_seq[which.min(tau_error)]
  lambda_min <- lambda_seq[which.min(lambda_error)]
  cat("best tau and lambda are : ", tau_min, lambda_min, "\n")
  # MLE result
  fit_ML <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = 0, standardize = F)
  beta_ML <- as.numeric(coef(fit_ML))

  # ridge result
  fit_ridge <- glmnet(X, Surv(TT, status), family = "cox", alpha = 0, lambda = lambda_min, standardize = F)
  beta_ridge <- as.numeric(coef(fit_ridge))
  # catlaytic result
  beta_syn <- mimic_glmnet_cox_MLE_syn(TT, status, X, my_syndata, tau_min)

  error_table <- c(ML.MSE = norm(beta_ML - beta_0, "2")^2 / p, ridge.MSE = norm(beta_ridge - beta_0, "2")^2 / p, syn.MSE = norm(beta_syn - beta_0, "2")^2 / p)
  list(error_table = error_table, tau_seq = tau_seq, tau_cv_error = tau_error, lambda_seq = lambda_seq, lambda_CV_error = lambda_error)
}


generate_beta_0 <- function(p = 20, r = 1, sparse_prop = 0.25) {
  rbinom(p, 1, prob = 1 - sparse_prop) * r
}
generate_covariate_cox <- function(n, p = 10) {
  tem_f <- function(p) {
    Zij_X <- rnorm(p)
    index_p <- (2 * (1:p)) < p
    X <- Zij_X
    X[index_p] <- 2 * (Zij_X[index_p] > 0) - 1

    return(X)
  }

  # tem_f=rnorm(p)
  design <- t(sapply(1:n, function(i) {
    tem_f(p)
  }))
  design
}
generate_covariate_cox2 <- function(n, p = 10) {
  design <- t(sapply(1:n, function(i) {
    rnorm(p)
  }))
  design
}




generate_synthetic_data_no_need_X <- function(M, TT, X, status, local_cnt = 5) {
  X.syn <- IRS_M_each_people_X(X, M)
  event_index <- which(status == 1)
  TT_event <- TT[event_index]
  T.syn <- sample(TT_event, M, replace = T)

  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}

generate_synthetic_data_pairdist_X <- function(M, TT, X, status, local_cnt = 5) {
  X.syn <- IRS_M_each_people_X(X, M)
  event_index <- which(status == 1)
  X_event <- X[event_index, ]
  TT_event <- TT[event_index]
  T.syn <- rep(0, M)
  # synthetic response is local average instead of global weighted average since some of value is big
  for (i in 1:M) {
    # whole_dist=sum(sapply(1:sum(status),function(j){norm(X_event[,j]-X.syn[i,],"2")}))
    pairdist <- sapply(1:sum(status), function(j) {
      norm(X_event[j, ] - X.syn[i, ], "2")
    })
    T.syn[i] <- mean(TT_event[order(pairdist)[1:local_cnt]])
  }
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}

generate_syn_response_based_on_history <- function(M, TT, X, status, beta_with_important_covariate) {
  X.syn <- IRS_M_each_people_X(X, M)
  event_index <- which(status == 1)
  X_event <- X[event_index, ]
  TT_event <- TT[event_index]
  T.syn <- rep(0, M)
  true_data_inner_prod_set <- as.numeric(X_event %*% beta_with_important_covariate) # length is sum(status)
  syn_data_inner_prod_set <- as.numeric(X.syn %*% beta_with_important_covariate) # length is M
  for (i in 1:M) {
    qi.star <- syn_data_inner_prod_set[i]
    dist <- sapply(1:length(true_data_inner_prod_set), function(index) {
      abs(qi.star - true_data_inner_prod_set[index])
    })
    order_dist <- order(dist) # now we can choose Top-20% close dist index, and then sample from these response uniformly
    top20_index <- order_dist[1:round(length(true_data_inner_prod_set) / 5)]
    # top20_index=order_dist[1]
    top20_candidate <- TT_event[top20_index]
    T.syn[i] <- sample(top20_candidate, 1)
  }
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}


generate_syn_response_based_on_history_percent <- function(M, TT, X, status, beta_with_important_covariate, top_percent) {
  X.syn <- IRS_M_each_people_X(X, M)
  event_index <- which(status == 1)
  X_event <- X[event_index, ]
  TT_event <- TT[event_index]
  T.syn <- rep(0, M)
  true_data_inner_prod_set <- as.numeric(X_event %*% beta_with_important_covariate) # length is sum(status)
  syn_data_inner_prod_set <- as.numeric(X.syn %*% beta_with_important_covariate) # length is M
  for (i in 1:M) {
    qi.star <- syn_data_inner_prod_set[i]
    dist <- sapply(1:length(true_data_inner_prod_set), function(index) {
      abs(qi.star - true_data_inner_prod_set[index])
    })
    order_dist <- order(dist) # now we can choose Top-20% close dist index, and then sample from these response uniformly
    top_index <- order_dist[1:round(length(true_data_inner_prod_set) * top_percent)]
    # top20_index=order_dist[1]
    top_candidate <- TT_event[top_index]
    T.syn[i] <- sample(top_candidate, 1)
  }
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}


generate_syn_response_based_on_history_spline <- function(M, TT, X, status, beta_with_important_covariate, top_percent) {
  X.syn <- IRS_M_each_people_X(X, M)
  event_index <- which(status == 1)
  X_event <- X[event_index, ]
  TT_event <- TT[event_index]
  T.syn <- rep(0, M)
  true_data_inner_prod_set <- as.numeric(X_event %*% beta_with_important_covariate) # length is sum(status)
  # apply smooth spline to fir y=g(x^T beta)
  sss <- smooth.spline(true_data_inner_prod_set, TT_event)
  syn_data_inner_prod_set <- as.numeric(X.syn %*% beta_with_important_covariate) # length is M
  T.syn <- predict(sss, syn_data_inner_prod_set)$y
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}

generate_synthetic_data_ridge_sub <- function(M, X, beta_sub) {
  X.syn <- IRS_M_each_people_X(X, M)
  T.syn <- sapply(1:M, function(i) {
    rexp(1, exp(as.numeric(X.syn %*% beta_sub))[i])
  })
  return(list(T.syn = T.syn, X.syn = X.syn, status.syn = rep(1, M)))
}

CoxGradient_dup <- function(beta, Xs_dup, entry, Ts_dup, event_dup, tau, M) {
  exp_Z_beta <- exp(Xs_dup %*% beta)
  n <- nrow(Xs_dup) - M
  weight <- c(rep(1, n), rep(tau / M, M))
  weighted_exp_Z_beta <- exp_Z_beta * weight
  p <- ncol(Xs_dup)

  grad_matrix_need_rowsum <- sapply(1:(M + n), function(i) {
    if (event_dup[i] == 0) {
      return(rep(0, p))
    }
    risk_index <- GetRiskSet(Ts_dup[i], entry, Ts_dup, event_dup)
    # print(risk_index)
    weighted_exp_Z_beta_risk <- weighted_exp_Z_beta[risk_index]
    # print(weighted_exp_Z_beta_risk)
    fenmu <- sum(weighted_exp_Z_beta_risk)
    if (length(risk_index) > 1) {
      X_risk_set <- Xs_dup[risk_index, ]
    } else {
      X_risk_set <- matrix(Xs_dup[risk_index, ], nc = ncol(Xs_dup))
    }
    if (i > n) {
      return(tau / M * (Xs_dup[i, ] - t(X_risk_set) %*% weighted_exp_Z_beta_risk / fenmu))
    } else {
      return(Xs_dup[i, ] - t(X_risk_set) %*% weighted_exp_Z_beta_risk / fenmu)
    }
  })
  # print(grad_matrix_need_rowsum)
  return(apply(grad_matrix_need_rowsum, 1, sum))
}

CoxHessian_dup <- function(beta, Xs_dup, entry, Ts_dup, event_dup, tau, M) { # slow is acceptable
  exp_Z_beta <- as.numeric(exp(Xs_dup %*% beta))
  n <- nrow(Xs_dup) - M
  weight <- c(rep(1, n), rep(tau / M, M))
  weighted_exp_Z_beta <- exp_Z_beta * weight
  p <- ncol(Xs_dup)


  matrix_list <- lapply(1:(M + n), function(i) {
    if (event_dup[i] == 0) {
      return(matrix(0, p, p))
    }
    risk_index <- GetRiskSet(Ts_dup[i], entry, Ts_dup, event_dup)
    weighted_exp_Z_beta_risk <- weighted_exp_Z_beta[risk_index]
    fenmu <- sum(weighted_exp_Z_beta_risk)
    if (length(risk_index) > 1) {
      X_risk_set <- Xs_dup[risk_index, ]
    } else {
      X_risk_set <- matrix(Xs_dup[risk_index, ], nc = ncol(Xs_dup))
    }
    tmp_sqrt_theta_X <- X_risk_set * sqrt(weighted_exp_Z_beta_risk)
    A1 <- t(tmp_sqrt_theta_X) %*% tmp_sqrt_theta_X * fenmu
    B0 <- t(X_risk_set) %*% weighted_exp_Z_beta_risk
    B1 <- B0 %*% t(B0)
    if (i > n) {
      return(tau / M * (A1 - B1) / fenmu^2)
    } else {
      return((A1 - B1) / fenmu^2)
    }
  })

  hess <- Reduce("+", matrix_list)

  return(-hess)
}



newton_cox_dup <- function(tol = 1e-5, MAX_ITER = 25, initial_values, TT, status, X, tau, my_syndata) {
  n <- length(TT)
  dup_X <- X
  dup_TT <- TT
  dup_status <- status
  dup_X <- rbind(dup_X, my_syndata$X.syn)
  dup_TT <- c(dup_TT, my_syndata$T.syn)
  dup_status <- c(dup_status, my_syndata$status.syn)
  M <- length(my_syndata$T.syn)

  entry <- rep(0, n + M)
  beta_old <- initial_values

  H_old <- CoxHessian_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
  gd_old <- CoxGradient_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
  beta_new <- beta_old - qr_solve_my(H_old, gd_old)
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    # print(beta_old)
    # print("here")
    beta_old <- beta_new
    H_old <- CoxHessian_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
    gd_old <- CoxGradient_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
    beta_new <- beta_old - qr_solve_my(H_old, gd_old)
    iter <- iter + 1
    # print(beta_new)
    if (is.na(beta_new[1])) {
      print("fail to compute, may try new initial value")
      return(initial_values + 10000)
    }
    if (norm(as.numeric(beta_new), "2") > 1000000) {
      print("does not converge")
      return(initial_values + 50000)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) cat("mixed not converge \n")
  # if(iter<MAX_ITER) cat("number of iteration to convergence: ",iter,"\n")
  return(as.numeric(beta_new))
}




qr_solve_my <- function(M, b) {
  # print(M)
  s <- qr.coef(qr(M), b)
  s[is.na(s)] <- 0.0001
  return(s)
}



newton_cox_syn <- function(tol = 1e-5, MAX_ITER = 25, initial_values, TT, status, X, tau, my_syndata) {
  n <- length(TT)
  entry <- rep(0, n)
  beta_old <- initial_values

  TT.star <- my_syndata$T.syn
  X.star <- my_syndata$X.syn
  status.star <- my_syndata$status.syn
  M <- length(TT.star)

  H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau / M * CoxHessian(beta_old, X.star, rep(0, M), TT.star, status.star)
  gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau / M * CoxGradient(beta_old, X.star, rep(0, M), TT.star, status.star)
  # print(H_old)
  # print(eigen(H_old)$value)
  # print(as.numeric(gd_old))
  beta_new <- beta_old - qr_solve_my(H_old, gd_old)
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    # print("here")
    beta_old <- beta_new
    H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau / M * CoxHessian(beta_old, X.star, rep(0, M), TT.star, status.star)
    # print(H_old[1:10,1:10])
    gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau / M * CoxGradient(beta_old, X.star, rep(0, M), TT.star, status.star)
    beta_new <- beta_old - qr_solve_my(H_old, gd_old)
    iter <- iter + 1
    # print(as.numeric(gd_old))
    # print(as.numeric(beta_new))

    if (is.na(beta_new[1])) {
      warning("NA BLEAK UP")
      return(initial_values + 10000)
    }
    if (norm(as.numeric(beta_new), "2") > 1000000) {
      warning("BLEAK UP")
      return(initial_values + 50000)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) warning("penalty not converge")
  # if(iter<MAX_ITER) cat("number of iteration to convergence: ",iter,"\n")
  return(as.numeric(beta_new))
}

newton_cox_ridge <- function(tol = 1e-5, MAX_ITER = 25, initial_values, TT, status, X, lambda) {
  n <- length(TT)
  entry <- rep(0, n)
  beta_old <- initial_values
  p <- length(beta_old)

  H_old <- CoxHessian(beta_old, X, entry, TT, status) - lambda * diag(p)
  # H_old=diag(diag(CoxHessian(beta_old,X,entry,TT,status)))+lambda*diag(p)

  gd_old <- CoxGradient(beta_old, X, entry, TT, status) - lambda * beta_old
  beta_new <- beta_old - solve(H_old, gd_old)
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    # print("here")
    beta_old <- beta_new
    H_old <- CoxHessian(beta_old, X, entry, TT, status) - lambda * diag(p)
    # H_old=diag(diag(CoxHessian(beta_old,X,entry,TT,status)))+lambda*diag(p)
    # print(diag(H_old))
    gd_old <- CoxGradient(beta_old, X, entry, TT, status) - lambda * beta_old
    beta_new <- beta_old - solve(H_old, gd_old)
    iter <- iter + 1
    # print(as.numeric(beta_new))

    if (is.na(beta_new[1])) {
      return(initial_values + 10000)
    }
    if (norm(as.numeric(beta_new), "2") > 1000000) {
      return(initial_values + 50000)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) cat("ridge not converge \n")
  # if(iter<MAX_ITER) cat("number of iteration to convergence: ",iter,"\n")
  return(as.numeric(beta_new))
}









CoxHessian_dup_onlydiag <- function(beta, Xs_dup, entry, Ts_dup, event_dup, tau, M) { # slow is acceptable
  exp_Z_beta <- as.numeric(exp(Xs_dup %*% beta))
  n <- nrow(Xs_dup) - M
  weight <- c(rep(1, n), rep(tau / M, M))
  weighted_exp_Z_beta <- exp_Z_beta * weight
  p <- ncol(Xs_dup)


  matrix_list <- lapply(1:(M + n), function(i) {
    if (event_dup[i] == 0) {
      return(rep(0, p))
    }
    risk_index <- GetRiskSet(Ts_dup[i], entry, Ts_dup, event_dup)
    weighted_exp_Z_beta_risk <- weighted_exp_Z_beta[risk_index]
    fenmu <- sum(weighted_exp_Z_beta_risk)
    if (length(risk_index) > 1) {
      X_risk_set <- Xs_dup[risk_index, ]
    } else {
      X_risk_set <- matrix(Xs_dup[risk_index, ], nc = ncol(Xs_dup))
    }
    tmp_sqrt_theta_X <- X_risk_set * sqrt(weighted_exp_Z_beta_risk)
    A1 <- t(tmp_sqrt_theta_X) %*% tmp_sqrt_theta_X * fenmu
    B0 <- t(X_risk_set) %*% weighted_exp_Z_beta_risk
    B1 <- B0 %*% t(B0)
    if (i > n) {
      return(diag(tau / M * (A1 - B1) / fenmu^2))
    } else {
      return(diag((A1 - B1) / fenmu^2))
    }
  })

  hess <- Reduce("+", matrix_list)

  return(-hess)
}


newton_cox_dup_hessonlydiag <- function(tol = 1e-5, MAX_ITER = 25, initial_values, TT, status, X, tau, my_syndata) {
  n <- length(TT)
  dup_X <- X
  dup_TT <- TT
  dup_status <- status
  dup_X <- rbind(dup_X, my_syndata$X.syn)
  dup_TT <- c(dup_TT, my_syndata$T.syn)
  dup_status <- c(dup_status, my_syndata$status.syn)
  M <- length(my_syndata$T.syn)

  entry <- rep(0, n + M)
  beta_old <- initial_values

  H_old <- CoxHessian_dup_onlydiag(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
  gd_old <- CoxGradient_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
  # beta_new=beta_old-qr_solve_my(H_old,gd_old)
  beta_new <- beta_old - gd_old / H_old
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    # print(beta_old)
    beta_old <- beta_new
    H_old <- CoxHessian_dup_onlydiag(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
    gd_old <- CoxGradient_dup(beta_old, dup_X, entry, dup_TT, dup_status, tau, M)
    beta_new <- beta_old - gd_old / H_old
    # beta_new=beta_old-qr_solve_my(H_old,gd_old)
    iter <- iter + 1
    # print(beta_new)
    if (is.na(beta_new[1])) {
      return(initial_values)
    }
    if (norm(as.numeric(beta_new), "2") > 10000) {
      return(initial_values)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) cat("not converge \n")
  # if(iter<MAX_ITER) cat("number of iteration to convergence: ",iter,"\n")
  return(as.numeric(beta_new))
}




newton_cox_syn_dec30_h_daga_const <- function(tol = 1e-5, MAX_ITER = 25, initial_values, TT, status, X, tau, my_syndata) {
  n <- length(TT)
  entry <- rep(0, n)
  beta_old <- initial_values

  TT.star <- my_syndata$T.syn
  X.star <- my_syndata$X.syn
  status.star <- my_syndata$status.syn
  h_daga <- my_syndata$h_star
  M <- length(TT.star)

  Xstarbetaold <- as.numeric(X.star %*% beta_old)
  H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau / M * CoxstarHessian(Xstarbetaold, X.star, TT.star, h_daga)
  gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau / M * CoxstarGradient(Xstarbetaold, X.star, TT.star, h_daga)
  # print(H_old)
  # print(eigen(H_old)$value)
  # print(as.numeric(gd_old))
  beta_new <- beta_old - qr_solve_my(H_old, gd_old)
  iter <- 1
  diff_new_old <- TRUE
  while ((iter < MAX_ITER) & diff_new_old) {
    # print("here")
    beta_old <- beta_new
    Xstarbetaold <- as.numeric(X.star %*% beta_old)
    H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau / M * CoxstarHessian(Xstarbetaold, X.star, TT.star, h_daga)
    # print(H_old[1:10,1:10])
    gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau / M * CoxstarGradient(Xstarbetaold, X.star, TT.star, h_daga)
    beta_new <- beta_old - qr_solve_my(H_old, gd_old)
    iter <- iter + 1
    # print(as.numeric(gd_old))
    # print(as.numeric(beta_new))

    if (is.na(beta_new[1])) {
      warning("NA BLEAK UP")
      return(initial_values + 10000)
    }
    if (norm(as.numeric(beta_new), "2") > 1000000) {
      warning("BLEAK UP")
      return(initial_values + 50000)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  if (iter == MAX_ITER) warning("penalty not converge")
  # if(iter<MAX_ITER) cat("number of iteration to convergence: ",iter,"\n")
  return(as.numeric(beta_new))
}

CoxstarGradient <- function(Xsbeta, Xs, Ts, h_daga) {
  # no need to add status since no censoring for synthetic data
  M <- nrow(Xs)
  p <- ncol(Xs)
  gd <- rep(0, p)
  for (ii in 1:M) {
    gd <- gd + (1 - Ts[ii] * h_daga * exp(Xsbeta[ii])) * Xs[ii, ]
  }
  return(gd)
}

CoxstarHessian <- function(Xsbeta, Xs, Ts, h_daga) {
  M <- nrow(Xs)
  p <- ncol(Xs)
  hess <- matrix(0, p, p)
  for (ii in 1:M) {
    hess <- hess - Ts[ii] * h_daga * exp(Xsbeta[ii]) * tcrossprod(Xs[ii, ])
  }
  return(hess)
}

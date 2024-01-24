frequentist_prior <-
  function(data_,
           taus = NULL,
           risk_estimate_method = "parametric_bootstrap",
           discrepancy_function = "square_error",
           bootstrap_iteration_times = 100,
           tau_0 = 1) {
    .not_in_stop(
      c("comb_y", "comb_x", "obs_n", "syn_n", "fit_fam", "obs_x", "obs_y", "reg_type", "x_col_n"),
      names(data_),
      "Not all necessary variables are present in data_."
    )

    if (is.null(taus)) taus <- data_$x_col_n / 4

    if (is.null(data_$est_beta)) data_$est_beta <- list()

    if (length(taus) == 1) {
      message("Will only get the point estimate for single tau value for faster process.")

      data_$est_beta <- unique(append(data_$est_beta, list(list(
        "method" = paste(
          "single point estimate:",
          "tau = ",
          taus
        ),
        "beta" = fit(data_$comb_x,
                     data_$comb_y,
                     weights = c(rep(1, data_$obs_n), rep(taus / data_$syn_n, data_$syn_n)),
                     family = data_$fit_fam
        )$coefficients))))

      return(data_)
    }

    if (risk_estimate_method == "mallows_estimate" &
        discrepancy_function != "square_error") stop("Please use square_error with mallows_estimate")

    min_risk_estimate <- Inf
    optimal_tau <- Inf
    for (tau in taus) {
      beta_tau <- fit(data_$comb_x,
        data_$comb_y,
        weights = c(rep(1, data_$obs_n), rep(tau / data_$syn_n, data_$syn_n)),
        family = data_$fit_fam
      )$coefficients

      D <- mean(discrepancy_functions(discrepancy_function)(data_$obs_x,
        data_$obs_y,
        beta_tau,
        data_$reg_type))

      omega <- do.call(
        risk_estimate_method,
        list(
          data_,
          tau,
          discrepancy_function,
          tau_0,
          bootstrap_iteration_times
        )
      )

      if (min_risk_estimate > D + omega) {
        min_risk_estimate <- D + omega
        optimal_tau <- tau
      }
    }
    data_$est_beta <- unique(append(data_$est_beta, list(list(
      "method" = paste(
        "frequentist:",
        "risk estimate method =",
        risk_estimate_method,
        "; discrepancy function =",
        discrepancy_function,
        "; bootstrap iteration times =",
        bootstrap_iteration_times,
        "; tau 0 = ",
        tau_0,
        "; optimal tau = ",
        optimal_tau
      ),
      "beta" = fit(data_$comb_x,
        data_$comb_y,
        weights = c(rep(1, data_$obs_n), rep(optimal_tau / data_$syn_n, data_$syn_n)),
        family = data_$fit_fam
      )$coefficients
    ))))


    return(data_)
  }


parametric_bootstrap <- function(data_,
                                 tau,
                                 discrepancy_function = "square_error",
                                 tau_0 = 1,
                                 bootstrap_iteration_times = 100) {
  .large_number_warning(data_$syn_n,
                        100,
                        "Number of synthetic data is large, suggest use cross-validation.")

  hat_tau_0_beta <- fit(data_$comb_x,
    data_$comb_y,
    weights = c(rep(1, data_$obs_n), rep(tau_0 / data_$syn_n, data_$syn_n)),
    family = data_$fit_fam
  )$coefficients

  mean_Y <- reg_mean(data_$reg_type)(data_$obs_x, hat_tau_0_beta)

  boot_Y_collect <-
    lambda_boot_hat_mu_collect <- matrix(NA,
      nrow = data_$obs_n,
      ncol = bootstrap_iteration_times
    )

  for (bootstrap_iter in 1:bootstrap_iteration_times) {
    boot_Y <- generate_data_from_dist(
      data_$fit_fam,
      data_$obs_n,
      mean_Y
    )

    boot_hat_tau_beta <- fit(data_$comb_x,
      c(boot_Y, data_$syn_y),
      weights = c(rep(1, data_$obs_n), rep(tau / data_$syn_n, data_$syn_n)),
      family = data_$fit_fam
    )$coefficients


    boot_Y_collect[, bootstrap_iter] <- boot_Y

    lambda_boot_hat_mu_collect[, bootstrap_iter] <- lambda_function(
      discrepancy_function
    )(data_$obs_x,
      boot_hat_tau_beta,
      data_$reg_type)
  }

  return(sum(sweep(boot_Y_collect, 1, rowMeans(boot_Y_collect)) *
    lambda_boot_hat_mu_collect) /
    bootstrap_iteration_times)
}

steinian_estimate <- function(data_,
                              tau,
                              discrepancy_function = "square_error",
                              tau_0 = 1,
                              ...) {
  .large_number_warning(data_$syn_n,
                        50,
                        "Number of synthetic data is large, suggest use cross-validation\\parametric_bootstrap.")

  .not_same_stop(data_$reg_type,
                 "logistic",
                 "Only logistic regression can use steinian_estimate.")

  hat_beta_0 <- fit(data_$comb_x,
    data_$comb_y,
    weights = c(rep(1, data_$obs_n), rep(tau_0 / data_$syn_n, data_$syn_n)),
    family = data_$fit_fam
  )$coefficients

  hat_mu_0 <- reg_mean(data_$reg_type)(data_$obs_x, hat_beta_0)

  hat_tau_beta <- fit(data_$comb_x,
    data_$comb_y,
    weights = c(rep(1, data_$obs_n), rep(tau / data_$syn_n, data_$syn_n)),
    family = data_$fit_fam
  )$coefficients

  steinian_sum <- 0
  for (obs in 1:data_$obs_n) {
    triangle_response <- data_$comb_y
    triangle_response[obs] <- 1 - data_$comb_y[obs]

    trigangle_hat_tau_beta <- fit(data_$comb_x,
      triangle_response,
      c(rep(1, data_$obs_n), rep(tau / data_$syn_n, data_$syn_n)),
      family = data_$fit_fam
    )$coefficients

    steinian_sum <- steinian_sum +
      hat_mu_0[obs] *
        (1 - hat_mu_0)[obs] *
        (2 * data_$comb_x[obs, 1] - 1) *
        (lambda_function(discrepancy_function)(
          data_$comb_x[obs, ], hat_tau_beta, data_$reg_type) -
          lambda_function(discrepancy_function)(
            data_$comb_x[obs, ], trigangle_hat_tau_beta, data_$reg_type)
        )
  }

  return(steinian_sum / data_$obs_n)
}

mallows_estimate <- function(data_,
                             tau,
                             discrepancy_function = "square_error",
                             tau_0 = 1,
                             ...) {

  .not_same_stop(data_$reg_type,
                 "linear",
                 "Only linear regression can use steinian_estimate.")

  beta_0 <- fit(data_$comb_x,
                data_$comb_y,
                weights = c(rep(1, data_$obs_n), rep(tau_0 / data_$syn_n, data_$syn_n)),
                family = data_$fit_fam
  )$coefficients

  var_Y <-

  tau_H <-

  return(2*var_Y*sum(diag(tau_H))/length(tau_H))
}

H_matrix <- function(X, X_star, tau, M, Y) {

  X <- as.matrix(X)
  X_star <- as.matrix(X_star)
  Y <- as.matrix(Y)
  return(
    t(X)%*%(X%*%t(X)+(tau/M)*t(X_star)%*%X_star)^(-1)%*%
      (t(X)+tau/(M*nrow(X))*t(X_star)%*%matrix(rep(1, M))%*%t(matrix(rep(1,nrow(X)))))
  )

}

# TODO: confirm it is belong to frequentist/ seperate
cross_validation <- function(data_,
                             tau,
                             discrepancy_function = "square_error",
                             tau_0 = 1,
                             cv_fold_num = 5,
                             ...
                             ) {

  error_sum <- 0

  for (fold in 1:cv_fold_num) {

    test_x <- data_$obs_x[(nrow(data_$obs_x)*(fold-1)/cv_fold_num+1):
                  (nrow(data_$obs_x)*fold/cv_fold_num),]
    test_y <- data_$obs_y[(nrow(data_$obs_x)*(fold-1)/cv_fold_num+1):
                            (nrow(data_$obs_x)*fold/cv_fold_num)]

    train_x <- data_$obs_x[-(nrow(data_$obs_x)*(fold-1)/cv_fold_num+1):
                                -(nrow(data_$obs_x)*fold/cv_fold_num),]
    train_y <- data_$obs_y[-(nrow(data_$obs_x)*(fold-1)/cv_fold_num+1):
                             -(nrow(data_$obs_x)*fold/cv_fold_num)]

    estimate_beta <- fit(train_x, train_y, family=data_$fit_fam)$coefficients

    error_sum <- error_sum + mean(discrepancy_functions(discrepancy_function)(test_x,
                                                                              test_y,
                                                                              estimate_beta,
                                                                              data_$reg_type))
  }

  return(error_sum/cv_fold_num)

}


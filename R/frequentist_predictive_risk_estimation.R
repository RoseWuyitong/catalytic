frequentist_predictive_risk_estimation <- # TODO: change name, e.g. tune tau
  function(data_,
           taus = NULL, # TODO: if user input one value, warning and select a range of taus around the tau.
           omega_estimate_method = "parametric_bootstrap",
           discrepancy_function = "square_error",
           bootstrap_iteration_times = 100,
           tau_0 = 1) {
    .not_in_stop(
      c("comb_y", "comb_x", "obs_n", "syn_n", "fit_fam", "obs_x", "obs_y", "reg_type"),
      names(data_),
      "Not all necessary variables are present in data_."
    )

    if (is.null(taus)) {
      taus <- data_$syn_n / 4 # should be p/4
    } else if (is.double(taus)) {
      taus <- taus
    } else {
      stop("Please enter taus in double format")
    }


    .large_number_warning(data_$syn_n,
                          100,
                          "Number of synthetic data is large, suggest use cross-validation.")

    if (is.null(data_$est_beta)) data_$est_beta <- list()

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
        omega_estimate_method,
        list(
          data_,
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
        "omega estimate method =",
        omega_estimate_method,
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
                                 discrepancy_function = "square_error",
                                 tau_0 = 1,
                                 bootstrap_iteration_times = 100) {
  # TODO: threhold should be 100
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
                              discrepancy_function = "square_error",
                              tau_0 = 1,
                              bootstrap_iteration_times = NULL) {
  # TODO: threhold should be 50
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

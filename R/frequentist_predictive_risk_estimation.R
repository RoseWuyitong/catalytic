frequentist_predictive_risk_estimation <-
  function(data,
           taus = 1:10,
           discrepancy_function = "square_error",
           omega_estimate_method = "parametric_bootstrap",
           boot_data_generation_dist = "uniform",
           boot_data_generation_dist_sd = 0,
           log_density_weights = FALSE,
           fit_family = "binomial",
           tau_0 = 1,
           bootstrap_iteration_times = 100) {
    min_risk_estimate <- Inf


    for (b in bootstrap_iteration_times) {
      for (tau in taus) {
        beta_tau <- fit(data$comb_x,
          data$comb_y,
          c(
            rep(1, data$obs_n),
            rep(tau / data$syn_n, data$syn_n)
          ),
          family = fit_family
        )$coefficients
      }
    }



    return(data)
  }

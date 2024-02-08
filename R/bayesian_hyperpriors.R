bayesian_hyperpriors <- function(data_,
                                 joint = FALSE,
                                 tau = NULL,
                                 chains = 4,
                                 iter = 2000,
                                 warmup = 1000,
                                 mcmc_algorithm = "HMC",
                                 joint_alpha = 2,
                                 joint_gamma = 1) {
  .not_in_stop(
    c("comb_y", "comb_x", "obs_n", "syn_n", "fit_fam", "obs_x", "obs_y", "reg_type", "x_col_n"),
    names(data_),
    "Not all necessary variables are present in data_."
  )

  if (is.null(tau) & !joint) tau <- data_$x_col_n / 4

  if (is.null(data_$est_beta)) data_$est_beta <- list()

  rstan_model <- rstan::stan(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      data_$reg_type,
      if (joint) "_joint",
      "_bayesian.stan"
    ),
    data = list(
      obs_n = data_$obs_n,
      syn_n = data_$syn_n,
      x_col_n = data_$x_col_n,
      obs_x = data_$obs_x,
      obs_y = data_$obs_y,
      syn_x = data_$syn_x,
      syn_y = data_$syn_y,
      tau = tau,
      kappa = log_density(data_$reg_type)(data_$syn_x,
                                          data_$syn_y,
                                          fit(data_$syn_x,
                                              data_$syn_y,
                                              family=data_$fit_fam)$coefficients,
                                          weights=1/data_$syn_n),
      alpha = joint_alpha,
      gamma = joint_gamma
    ),
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = mcmc_algorithm
  )
  data_$est_beta <- unique(append(data_$est_beta, list(c(
    list(
      "method" = paste(
        "bayesian:",
        "joint =",
        joint,
        "; tau = ",
        tau,
        "; chains = ",
        chains,
        "; iter = ",
        iter,
        "; warmup = ",
        warmup
      ),
      "rhat_pass" = gelman_rubin_statistic(rstan_model),
      "ess_pass" = effective_sample_size(rstan_model),
      "rhat" = head(rstan::summary(rstan_model)$summary[, "Rhat"], -1),
      "ess" = head(rstan::summary(rstan_model)$summary[, "n_eff"], -1),
      "model" = rstan_model
    ),
    if (joint) {
      list(
        "beta_median" = as.numeric(head(rstan::summary(rstan_model)$summary[, "50%"][-1], -1)),
        "beta_lower_iqr" = as.numeric(head(rstan::summary(rstan_model)$summary[, "2.5%"][-1], -1)),
        "beta_higher_iqr" = as.numeric(head(rstan::summary(rstan_model)$summary[, "97.5%"][-1], -1)),
        "beta_mean" = as.numeric(head(rstan::summary(rstan_model)$summary[, "mean"][-1], -1)),
        "tau_median" = as.numeric(rstan::summary(rstan_model)$summary[, "50%"][1]),
        "tau_lower_iqr" = as.numeric(rstan::summary(rstan_model)$summary[, "2.5%"][1]),
        "tau_higher_iqr" = as.numeric(rstan::summary(rstan_model)$summary[, "97.5%"][1]),
        "tau_mean" = as.numeric(rstan::summary(rstan_model)$summary[, "mean"][1])
      )
    } else {
      list(
        "beta_median" = as.numeric(head(rstan::summary(rstan_model)$summary[, "50%"], -1)),
        "beta_lower_iqr" = as.numeric(head(rstan::summary(rstan_model)$summary[, "2.5%"], -1)),
        "beta_higher_iqr" = as.numeric(head(rstan::summary(rstan_model)$summary[, "97.5%"], -1)),
        "beta_mean" = as.numeric(head(rstan::summary(rstan_model)$summary[, "mean"], -1))
      )
    }
  ))))

  return(data_)
}

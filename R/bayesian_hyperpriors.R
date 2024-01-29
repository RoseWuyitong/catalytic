bayesian_hyperpriors <- function(data_,
                                 taus = NULL,
                                 joint = FALSE) {
  .not_in_stop(
    c("comb_y", "comb_x", "obs_n", "syn_n", "fit_fam", "obs_x", "obs_y", "reg_type", "x_col_n"),
    names(data_),
    "Not all necessary variables are present in data_."
  )

  if (!joint & is.null(tau)) tau <- data_$x_col_n / 4

  data_ <- ifelse(joint, joint_bayesian(data_, taus), not_joint_bayesian(data_, taus))
}


joint_bayesian <- function() {

}

# Function to perform Gibbs sampling
gibbs_sampling <- function(data_,
                           taus = NULL,
                           tau_0 = 1,
                           alpha = 2,
                           gamma = 1,
                           iteration_times = 100) {
  if (is.null(taus)) taus <- seq(max(data_$x_col_n / 4 - 2, 1),
                                 data_$x_col_n / 4 + 2,
                                 by = 0.5)


  gibbs_samples <- matrix(NA, nrow = iteration_times , ncol = data_$x_col_n + 2)
  colnames(samples) <- c("tau", paste0("beta_", 0:data_$x_col_n))

  gibbs_samples[1, ] <- c(tau_0, fit(data_$comb_x,
    data_$comb_y,
    weights = c(
      rep(1, data_$obs_n),
      rep(tau_0, data_$syn_n)
    ),
    family = data_$fit_fam
  )$coefficients)

  for (iter in 2:iteration_times) {
    # update tau
    kappa <- log_density(data_$reg_type)(data_$syn_x,
                                             data_$syn_y,
                                             fit(data_$syn_x,
                                                 data_$syn_y,
                                                 family=data_$fit_fam)$coefficients,
                                         weights=1/data_$syn_n)

    c_value <- ifelse(identical(data_$reg_type, "linear"), (data_$x_col_n + alpha) / 2, data_$x_col_n + alpha)

    taus_prob <- taus^(c_value-1) * exp(tau*(
      log_density(data_$reg_type)(data_$syn_x,
      data_$syn_y,
      gibbs_samples[iter-1, -1], weights=1/data_$syn_n) - (
        kappa + 1/gamma)))

    tau <- taus[sample(1:length(taus), 1, prob = taus_prob)]

    # update beta by HMC
    beta <- hamiltonian_monte_carlo(X, Y, gibbs_samples[iter-1, -1], data_$reg)

    # Save samples
    gibbs_samples[iter, ] <- c(tau, beta)
  }

  return(samples)
}


hamiltonian_monte_carlo <- function(
    X,
    Y,
    beta,
    reg_type,
    adptive_variance = 1
    ) {
  step_size <- rnorm(1)

  adaptive_variance <- (pracma::hessian(-log_density(reg_type)(X, Y, beta)))^(-1)
}
















not_joint_bayesian <- function(data_,
                               tau= NULL,
                               tau_0 = 1,
                               iteration_times = 100,
                               ...) {

  beta_0 <- fit(data_$comb_x,
                data_$comb_y,
                weights = c(
                  rep(1, data_$obs_n),
                  rep(tau_0, data_$syn_n)
                ),
                family = data_$fit_fam
  )$coefficients

  post_density <- exp(log_density(data_$reg_type)(data_$syn_x,
                                                  data_$syn_y,
                                                  beta_0,
                                                  weights=tau/data_$syn_n)+
                        log_density(data_$reg_type)(data_$obs_x,
                                                    data_$obs_y,
                                                    beta_0))

  for (iter in 2:iteration_times) {}

}


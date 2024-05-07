start_a <- function(...) {
  Z_corr <- 0.5
  Z_mean <- 0
  Z_variance <- 1

  z <- 0.75 # nonsparsity
  r <- 0.3 # oracle classification error

  sample_size <- 30 # observed data sample size CHANGE
  dimention <- 20 # number of covariate CHNAGE
  M <- 200

  init_taus <- seq(1, 100, 2) # CHANGE
  simulation_times <- 3 # CHANGE

  iter <- 7000
  warmup <- 3000

  norm_stan_model_file <- "~/catalytic/R/rstan_script/binomial_joint_bayesian_theta_norm.stan"
  square_norm_stan_model_file <- "~/catalytic/R/rstan_script/binomial_joint_bayesian_theta_norm_square.stan"

  gamma_stan_model_file <- "~/catalytic/R/rstan_script/binomial_joint_bayesian_theta_gamma.stan"
  multi_gamma_stan_model_file <- "~/catalytic/R/rstan_script/binomial_joint_bayesian_theta_gamma_multi.stan"


  glm_results <-
    freq_boots_results <-
    freq_stein_results <-
    freq_cv_results <-
    bayes_tau_01_nuts_mean_results <-
    bayes_theta_2_nuts_mean_results <-
    bayes_theta_2pn_nuts_mean_results <-
    bayes_theta_joint_alpha_norm_nuts_mean_results <-
    bayes_theta_joint_alpha_norm_square_nuts_mean_results <-
    bayes_theta_joint_alpha_gamma_nuts_mean_results <-
    bayes_theta_joint_alpha_gamma_multi_nuts_mean_results <-
    bayes_gibbs_5_mean_results <-
    bayes_gibbs_5_joint_alpha_mean_results <-
    bayes_gibbs_5_boots_mean_results <-
    glm_times <-
    freq_boots_times <-
    freq_stein_times <-
    freq_cv_times <-
    bayes_tau_01_nuts_times <-
    bayes_theta_2_nuts_times <-
    bayes_theta_2pn_nuts_times <-
    bayes_theta_joint_alpha_norm_nuts_times <-
    bayes_theta_joint_alpha_norm_square_nuts_times <-
    bayes_theta_joint_alpha_gamma_nuts_times <-
    bayes_theta_joint_alpha_gamma_multi_nuts_times <-
    bayes_gibbs_5_times <-
    bayes_gibbs_5_joint_alpha_times <-
    bayes_gibbs_5_boots_times <-
    c()

  suppressWarnings(
    for (t in (1:simulation_times)) {
      set.seed(t)
      cat("\nsimulation times:", t)

      observed_covariate <- observed_covariate_vector(
        dimention,
        sample_size,
        Z_mean,
        diag(rep(Z_variance, dimention)),
        Z_corr
      )

      real_coefficients <- real_coefficients_vector(dimention, z, r)

      observed_response <- data.frame(y = get_data("binomial")(
        n = sample_size,
        mean = get_mean("binomial")(X = observed_covariate,
          coefficients = real_coefficients)
      ))
      obs_data <- cbind(observed_covariate, observed_response)

      data_ <- data_initialization(y ~ 1, family = binomial(), data = obs_data, syn_size = M)


      cat("\nglm")
      start_time <- Sys.time()
      glm_coefficients <- stats::glm(y ~ ., family = binomial(), data = obs_data)$coefficients
      glm_time <- difftime(Sys.time(), start_time, units = "secs")

      cat("\nfreq_boots")
      start_time <- Sys.time()
      freq_boots <- cat_glm_tune(y ~ .,
        family = binomial(),
        data_ = data_,
        risk_estimate_method = "parametric_bootstrap",
        discrepancy_function = "logistic_deviance",
        taus = init_taus
      )
      freq_boots_coefficients <- freq_boots$coefficients
      freq_boots_time <- difftime(Sys.time(), start_time, units = "secs")

      cat("\nfreq_stein_coefficients")
      start_time <- Sys.time()
      freq_stein <- cat_glm_tune(y ~ .,
        family = binomial(),
        data_ = data_,
        risk_estimate_method = "steinian_estimate",
        discrepancy_function = "logistic_deviance",
        taus = init_taus
      )
      freq_stein_coefficients <- freq_stein$coefficients
      freq_stein_time <- difftime(Sys.time(), start_time, units = "secs")

      cat("\nfreq_cv")
      start_time <- Sys.time()
      freq_cv <- cat_glm_tune(y ~ .,
        family = binomial(),
        data_ = data_,
        risk_estimate_method = "cross_validation",
        discrepancy_function = "logistic_deviance",
        taus = init_taus
      )
      freq_cv_coefficients <- freq_cv$coefficients
      freq_cv_time <- difftime(Sys.time(), start_time, units = "secs")


      cat("\nbayes_tau_01_NUTS_model")
      start_time <- Sys.time()
      bayes_tau_01_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        iter = iter,
        warmup = warmup,
        tau_lower = 0.1,
        algorithm = "NUTS"
      )
      bayes_tau_01_nuts_time <- difftime(Sys.time(), start_time, units = "secs")

      bayes_tau_01_nuts_mean_coefficients <- bayes_tau_01_nuts_model$coefficients

      cat("\nbayes_theta_2_nuts_model")
      start_time <- Sys.time()
      bayes_theta_2_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        alpha = 2,
        algorithm = "NUTS"
      )
      bayes_theta_2_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_2_nuts_mean_coefficients <- bayes_theta_2_nuts_model$coefficients



      cat("\nbayes_theta_2pn_nuts_model")
      start_time <- Sys.time()
      bayes_theta_2pn_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        alpha = 2 + 2 * dimention / sample_size,
        algorithm = "NUTS"
      )
      bayes_theta_2pn_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_2pn_nuts_mean_coefficients <- bayes_theta_2pn_nuts_model$coefficients



      cat("\nbayes_theta_joint_alpha_norm_nuts_model")
      start_time <- Sys.time()
      bayes_theta_joint_alpha_norm_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        stan_model_file = norm_stan_model_file,
        algorithm = "NUTS"
      )
      bayes_theta_joint_alpha_norm_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_joint_alpha_norm_nuts_mean_coefficients <- bayes_theta_joint_alpha_norm_nuts_model$coefficients

      cat("\nbayes_theta_joint_alpha_gamma_nuts_model")
      start_time <- Sys.time()
      bayes_theta_joint_alpha_gamma_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        stan_model_file = gamma_stan_model_file,
        algorithm = "NUTS"
      )
      bayes_theta_joint_alpha_gamma_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_joint_alpha_gamma_nuts_mean_coefficients <- bayes_theta_joint_alpha_gamma_nuts_model$coefficients


      cat("\nbayes_theta_joint_alpha_norm_square_nuts_model")
      start_time <- Sys.time()
      bayes_theta_joint_alpha_norm_square_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        stan_model_file = square_norm_stan_model_file,
        algorithm = "NUTS"
      )
      bayes_theta_joint_alpha_norm_square_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_joint_alpha_norm_square_nuts_mean_coefficients <- bayes_theta_joint_alpha_norm_square_nuts_model$coefficients

      cat("\nbayes_theta_joint_alpha_gamma_multi_nuts_model")
      start_time <- Sys.time()
      bayes_theta_joint_alpha_gamma_multi_nuts_model <- cat_glm_bayes_joint(
        y ~ .,
        family = binomial,
        data_ = data_,
        use_theta = TRUE,
        iter = iter,
        warmup = warmup,
        stan_model_file = multi_gamma_stan_model_file,
        algorithm = "NUTS"
      )
      bayes_theta_joint_alpha_gamma_multi_nuts_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_theta_joint_alpha_gamma_multi_nuts_mean_coefficients <- bayes_theta_joint_alpha_gamma_multi_nuts_model$coefficients


      cat("cat_joint_bayesian_gibbs_glm")
      start_time <- Sys.time()
      bayes_gibbs_5_model <- cat_glm_bayes_joint_gibbs(
        y ~ .,
        family = binomial,
        data_ = data_,
        gibbs_iter = iter,
        gibbs_warmup = warmup
      )
      bayes_gibbs_5_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_gibbs_5_mean_coefficients <- bayes_gibbs_5_model$inform_df[-1, "mean"]


      cat("cat_joint_bayesian_gibbs_joint_alpha_glm")
      start_time <- Sys.time()
      bayes_gibbs_5_joint_alpha_model <- cat_glm_bayes_joint_gibbs(
        y ~ .,
        family = binomial,
        data_ = data_,
        gibbs_iter = iter,
        gibbs_warmup = warmup,
        init_tau = bayes_theta_joint_alpha_norm_nuts_model$tau,
        init_coefficients = bayes_theta_joint_alpha_norm_nuts_model$coefficients
      )
      bayes_gibbs_5_joint_alpha_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_gibbs_5_joint_alpha_mean_coefficients <- bayes_gibbs_5_joint_alpha_model$inform_df[-1, "mean"]

      cat("cat_joint_bayesian_gibbs_boots_glm")
      start_time <- Sys.time()
      bayes_gibbs_5_boots_model <- cat_glm_bayes_joint_gibbs(
        y ~ .,
        family = binomial,
        data_ = data_,
        gibbs_iter = iter,
        gibbs_warmup = warmup,
        init_tau = freq_boots$tau,
        init_coefficients = freq_boots$coefficients
      )
      bayes_gibbs_5_boots_time <- difftime(Sys.time(), start_time, units = "secs")
      bayes_gibbs_5_boots_mean_coefficients <- bayes_gibbs_5_boots_model$inform_df[-1, "mean"]


      test_x <- observed_covariate_vector(
        dimention,
        1000,
        Z_mean,
        diag(rep(Z_variance, dimention)),
        Z_corr
      )

      glm_results <- append(
        glm_results,
        err(test_x, real_coefficients, glm_coefficients)
      )
      freq_boots_results <- append(
        freq_boots_results,
        err(test_x, real_coefficients, freq_boots_coefficients)
      )
      freq_stein_results <- append(
        freq_stein_results,
        err(test_x, real_coefficients, freq_stein_coefficients)
      )
      freq_cv_results <- append(
        freq_cv_results,
        err(test_x, real_coefficients, freq_cv_coefficients)
      )

      bayes_tau_01_nuts_mean_results <- append(
        bayes_tau_01_nuts_mean_results,
        err(test_x, real_coefficients, bayes_tau_01_nuts_mean_coefficients)
      )

      bayes_theta_2_nuts_mean_results <- append(
        bayes_theta_2_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_2_nuts_mean_coefficients)
      )
      bayes_theta_2pn_nuts_mean_results <- append(
        bayes_theta_2pn_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_2pn_nuts_mean_coefficients)
      )

      bayes_theta_joint_alpha_norm_nuts_mean_results <- append(
        bayes_theta_joint_alpha_norm_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_joint_alpha_norm_nuts_mean_coefficients)
      )
      bayes_theta_joint_alpha_norm_square_nuts_mean_results <- append(
        bayes_theta_joint_alpha_norm_square_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_joint_alpha_norm_square_nuts_mean_coefficients)
      )

      bayes_theta_joint_alpha_gamma_nuts_mean_results <- append(
        bayes_theta_joint_alpha_gamma_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_joint_alpha_gamma_nuts_mean_coefficients)
      )
      bayes_theta_joint_alpha_gamma_multi_nuts_mean_results <- append(
        bayes_theta_joint_alpha_gamma_multi_nuts_mean_results,
        err(test_x, real_coefficients, bayes_theta_joint_alpha_gamma_multi_nuts_mean_coefficients)
      )


      bayes_gibbs_5_mean_results <- append(
        bayes_gibbs_5_mean_results,
        err(test_x, real_coefficients, bayes_gibbs_5_mean_coefficients)
      )

      bayes_gibbs_5_joint_alpha_mean_results <- append(
        bayes_gibbs_5_joint_alpha_mean_results,
        err(test_x, real_coefficients, bayes_gibbs_5_joint_alpha_mean_coefficients)
      )

      bayes_gibbs_5_boots_mean_results <- append(
        bayes_gibbs_5_boots_mean_results,
        err(test_x, real_coefficients, bayes_gibbs_5_boots_mean_coefficients)
      )

      glm_times <- append(
        glm_times,
        glm_time
      )
      freq_boots_times <- append(
        freq_boots_times,
        freq_boots_time
      )
      freq_stein_times <- append(
        freq_stein_times,
        freq_stein_time
      )
      freq_cv_times <- append(
        freq_cv_times,
        freq_cv_time
      )

      bayes_tau_01_nuts_times <- append(
        bayes_tau_01_nuts_times,
        bayes_tau_01_nuts_time
      )

      bayes_theta_2_nuts_times <- append(
        bayes_theta_2_nuts_times,
        bayes_theta_2_nuts_time
      )

      bayes_theta_2pn_nuts_times <- append(
        bayes_theta_2pn_nuts_times,
        bayes_theta_2pn_nuts_time
      )

      bayes_theta_joint_alpha_norm_nuts_times <- append(
        bayes_theta_joint_alpha_norm_nuts_times,
        bayes_theta_joint_alpha_norm_nuts_time
      )
      bayes_theta_joint_alpha_norm_square_nuts_times <- append(
        bayes_theta_joint_alpha_norm_square_nuts_times,
        bayes_theta_joint_alpha_norm_square_nuts_time
      )

      bayes_theta_joint_alpha_gamma_nuts_times <- append(
        bayes_theta_joint_alpha_gamma_nuts_times,
        bayes_theta_joint_alpha_gamma_nuts_time
      )
      bayes_theta_joint_alpha_gamma_multi_nuts_times <- append(
        bayes_theta_joint_alpha_gamma_multi_nuts_times,
        bayes_theta_joint_alpha_gamma_multi_nuts_time
      )

      bayes_gibbs_5_times <- append(
        bayes_gibbs_5_times,
        bayes_gibbs_5_time
      )

      bayes_gibbs_5_joint_alpha_times <- append(
        bayes_gibbs_5_joint_alpha_times,
        bayes_gibbs_5_joint_alpha_time
      )

      bayes_gibbs_5_boots_times <- append(
        bayes_gibbs_5_boots_times,
        bayes_gibbs_5_boots_time
      )
    }
  )


  # Define the list of methods
  methods <- c(
    "glm",
    "freq_boots",
    "freq_stein",
    "freq_cv",
    "bayes_tau_01",
    "bayes_theta_2",
    "bayes_theta_2pn",
    "bayes_theta_joint_alpha_norm",
    "bayes_theta_joint_alpha_norm_square",
    "bayes_theta_joint_alpha_gamma",
    "bayes_theta_joint_alpha_gamma_multi",
    "bayes_gibbs_5",
    "bayes_gibbs_joint_alpha_5",
    "bayes_gibbs_boots_5"
  )

  times_list <- list(
    glm_times,
    freq_boots_times,
    freq_stein_times,
    freq_cv_times,
    bayes_tau_01_nuts_times,
    bayes_theta_2_nuts_times,
    bayes_theta_2pn_nuts_times,
    bayes_theta_joint_alpha_norm_nuts_times,
    bayes_theta_joint_alpha_norm_square_nuts_times,
    bayes_theta_joint_alpha_gamma_nuts_times,
    bayes_theta_joint_alpha_gamma_multi_nuts_times,
    bayes_gibbs_5_times,
    bayes_gibbs_5_joint_alpha_times,
    bayes_gibbs_5_boots_times
  )

  # Define the list of results
  results_list <- list(
    glm_results,
    freq_boots_results,
    freq_stein_results,
    freq_cv_results,
    bayes_tau_01_nuts_mean_results,
    bayes_theta_2_nuts_mean_results,
    bayes_theta_2pn_nuts_mean_results,
    bayes_theta_joint_alpha_norm_nuts_mean_results,
    bayes_theta_joint_alpha_norm_square_nuts_mean_results,
    bayes_theta_joint_alpha_gamma_nuts_mean_results,
    bayes_theta_joint_alpha_gamma_multi_nuts_mean_results,
    bayes_gibbs_5_mean_results,
    bayes_gibbs_5_joint_alpha_mean_results,
    bayes_gibbs_5_boots_mean_results
  )

  # Calculate mean and standard error for each method
  means <- sapply(results_list, mean)
  ses <- sapply(results_list, function(x) sd(x) / sqrt(length(x)))
  mean_time <- round(sapply(times_list, mean), 2)
  ses_time <- round(sapply(times_list, function(x) sd(x) / sqrt(length(x))), 2)

  # Create the dataframe
  results <- data.frame(
    Method = methods,
    Mean = means,
    SE = ses,
    mean_time = mean_time,
    ses_time = ses_time
  )
}

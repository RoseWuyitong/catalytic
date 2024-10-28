#' Fit Bayesian Joint Cox Model with `cat_cox_bayes_joint`
#'
#' This function fits a Bayesian joint Cox proportional hazards model using Stan.
#' The model jointly estimates the baseline hazard and regression coefficients.
#'
#' @param formula A formula specifying the survival model. Must use `Surv` for survival analysis.
#' @param cat_init An object generated from `catalytic::cat_cox_initialization` containing initialization data.
#' @param ... Additional arguments passed to the `rstan::sampling` function.
#' @param stan_chains Number of Markov Chains (default is 4).
#' @param stan_iter Total number of iterations (default is 2000).
#' @param stan_warmup Number of warmup iterations (default is 1000).
#' @param tau_alpha Hyperparameter for the Gamma prior (default is 2).
#' @param tau_gamma Hyperparameter for the Gamma prior (default is 1).
#' @param hazard_beta Small value of `hazard_beta` corresponds to a more diffuse and less informative prior (default is 2).
#' @param H_star_list Optional list of cumulative hazard values. If NULL, it is estimated from a Weibull regression model.
#'
#' @return A list containing:
#' \item{function_name}{The name of the function.}
#' \item{formula}{The formula used for fitting the model.}
#' \item{cat_init}{The `cat_init` object used for initialization.}
#' \item{stan_chains}{Number of chains used in Stan sampling.}
#' \item{stan_iter}{Number of iterations used in Stan sampling.}
#' \item{stan_warmup}{Number of warmup iterations used in Stan sampling.}
#' \item{hazard_beta}{The prior parameter value `hazard_beta`.}
#' \item{H_star_list}{List of cumulative hazard values.}
#' \item{stan_model}{The compiled Stan model object.}
#' \item{stan_sample_model}{The fitted Stan model object.}
#' \item{coefficients}{Estimated regression coefficients.}
#' \item{increment_cumulative_baseline_hazard}{Estimated cumulative baseline hazard increments.}
#'
#' @export
#'
#' @examples
cat_cox_bayes_joint <- function(
    formula,
    cat_init,
    hazard_beta = 2,
    tau_alpha = 2,
    tau_gamma = 1,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    ...) {
  # Validate Input Parameters
  validate_cox_input(
    formula = formula,
    cat_init = cat_init,
    hazard_beta = hazard_beta,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  # Get unique times from sorted observed status times and a value larger than the maximum observed time
  unique_observed_status_times_and_beyond <- unique(c(
    # Sort the times where event occurred, status = 1
    sort(cat_init$obs_time[cat_init$obs_status == 1]),
    # Create a value larger than the maximum observed time
    2 * max(cat_init$obs_time) -
      max(cat_init$obs_time[-which(cat_init$obs_time == max(cat_init$obs_time))])
  ))


  # get risk set and failure set
  risk_and_failure_sets <- get_cox_risk_and_failure_sets(
    time_vector = cat_init$obs_time,
    status_vector = cat_init$obs_status,
    indicator_vector = unique_observed_status_times_and_beyond
  )


  # Get the initial guess of H_0 from the intercept-only weibull regression model based on observation data
  weibull_model <- do.call(
    survival::survreg,
    list(
      formula = survival::Surv(cat_init$obs_time, cat_init$obs_status) ~ 1,
      dist = "weibull"
    )
  )
  # Extracts the shape and scale parameters of the Weibull distribution
  weibull_shape <- as.numeric(1 / weibull_model$scale)
  weibull_scale <- as.numeric(exp(coef(weibull_model)))

  # Calculates the cumulative hazard values for the risk_and_failure_sets
  H_star_list <- (unique_observed_status_times_and_beyond / weibull_scale)^(weibull_shape)


  # List of tau_alpha value from each h_j ~ Gamma(alpha_j, hazard_beta), where alpha_j = alpha_0[j] - alpha_0[j-1]
  hazard_alpha_list <- diff(c(0, hazard_beta * H_star_list))

  stan_data <- list(
    obs_size = cat_init$obs_size,
    syn_size = cat_init$syn_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_x = cat_init$adj_syn_x,
    syn_time = cat_init$syn_time,
    syn_status = cat_init$syn_status,
    time_interval_num = length(unique_observed_status_times_and_beyond),
    adj_risk_set = risk_and_failure_sets$risk_set - risk_and_failure_sets$failure_set,
    failure_set = risk_and_failure_sets$failure_set,
    hazard_constant = cat_init$hazard_constant,
    hazard_alpha_list = hazard_alpha_list,
    hazard_beta = hazard_beta,
    kappa = get_cox_kappa(
      X = cat_init$syn_x,
      time = cat_init$syn_time,
      status = cat_init$syn_status,
      hazard_constant = cat_init$hazard_constant,
      entry_points = cat_init$entry_points
    ),
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  stan_model <- get_stan_model(
    type = "cox",
    joint_tau = TRUE
  )

  # Perform Bayesian Sampling via RStan
  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )

  ## Extract the mean of coefs from the stan sampling model
  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]

  if (!is.null(coefs)) {
    names(coefs) <- colnames(cat_init$adj_x)
  }

  cat_model <- list(
    function_name = "cat_cox_bayes_joint",
    ## Input/Processed parameters
    formula = as.formula(paste0(
      "survival::Surv(",
      cat_init$time_col_name,
      ",",
      cat_init$status_col_name,
      ")",
      f_rhs
    )),
    cat_init = cat_init,
    hazard_beta = hazard_beta,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...,
    ## Result
    model = stan_model,
    sample_model = stan_sample_model,
    tau = rstan::summary(stan_sample_model)$summary[
      grep("tau", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ],
    coefficients = coefs,
    increment_cumulative_baseline_hazard = rstan::summary(stan_sample_model)$summary[
      grep("h_j", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_model) <- c(cat_model$class, "cat_cox_bayes_joint", "bayes")


  return(cat_model)
}

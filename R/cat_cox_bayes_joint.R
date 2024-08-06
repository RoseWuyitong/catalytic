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
#' @param alpha Hyperparameter for the Gamma prior (default is 2).
#' @param gamma Hyperparameter for the Gamma prior (default is 1).
#' @param c_0 Small value of `c_0` corresponds to a more diffuse and less informative prior (default is 2).
#' @param H_star_list Optional list of cumulative hazard values. If NULL, it is estimated from a Weibull regression model.
#'
#' @return A list containing:
#' \item{function_name}{The name of the function.}
#' \item{formula}{The formula used for fitting the model.}
#' \item{cat_init}{The `cat_init` object used for initialization.}
#' \item{stan_chains}{Number of chains used in Stan sampling.}
#' \item{stan_iter}{Number of iterations used in Stan sampling.}
#' \item{stan_warmup}{Number of warmup iterations used in Stan sampling.}
#' \item{c_0}{The prior parameter value `c_0`.}
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
    ...,
    stan_chains = 4,
    stan_iter = 2000,
    stan_warmup = 1000,
    alpha = 2,
    gamma = 1,
    c_0 = 2,
    H_star_list = NULL) {
  if (!is.cat_cox_init(cat_init)) {
    stop(
      "Please enter `cat_init` that generated from `catalytic::cat_cox_initialization`."
    )
  }

  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  # Extract time, status and other columns(x) from `data`
  f_lhs <- get_formula_lhs(formula)
  if (f_lhs != "NULL" & !is.formula_surv(formula)) {
    stop("Please use `Surv` for formula.")

    f_lhs <- get_formula_lhs(formula)

    time_col_name <- gsub(
      "`",
      "",
      regmatches(
        f_lhs,
        regexec(
          "\\(([^,]+),",
          f_lhs
        )
      )[[1]][2]
    )
    status_col_name <- gsub(
      "`",
      "",
      regmatches(
        f_lhs,
        regexec(
          ",\\unique_observed_status_times_and_beyond*([^\\)]+)\\)",
          f_lhs
        )
      )[[1]][2]
    )

    if (time_col_name == cat_init$time_col_name &
      status_col_name == cat_init$status_col_name) {
      stop("Please use same time and status column name for formula")
    }
  }

  cat_init <- get_adjusted_data(cat_init = cat_init, formula = formula)

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
    time_vector = cat_init$obs_time[, 1],
    status_vector = cat_init$obs_status[, 1],
    indicator_vector = unique_observed_status_times_and_beyond
  )

  if (is.null(H_star_list)) {
    # Get the initial guess of H_0 from the intercept-only weibull regression model based on observation data
    surv_weibull_model <- do.call(
      survival::survreg,
      list(
        formula = survival::Surv(cat_init$obs_time[, 1], cat_init$obs_status[, 1]) ~ 1,
        dist = "weibull"
      )
    )
    # Extracts the shape and scale parameters of the Weibull distribution
    weibull_shape <- as.numeric(1 / surv_weibull_model$scale)
    weibull_scale <- as.numeric(exp(coef(surv_weibull_model)))

    # Calculates the cumulative hazard values for the risk_and_failure_sets
    H_star_list <- (unique_observed_status_times_and_beyond / weibull_scale)^(weibull_shape)
  }


  # List of alpha_0 value
  alpha_0_list <- c_0 * H_star_list

  # List of alpha value from each h_j ~ Gamma(alpha_j, c_0), where alpha_j = alpha_0[j] - alpha_0[j-1]
  alpha_of_h_j_list <- diff(c(0, alpha_0_list))

  stan_model <- get_stan_model(
    type = "cox",
    joint_tau = TRUE
  )

  stan_data <- list(
    J = length(unique_observed_status_times_and_beyond),
    alpha_of_h_j_list = alpha_of_h_j_list,
    c_0 = c_0,
    obs_size = cat_init$obs_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_size = cat_init$syn_size,
    syn_x = cat_init$adj_syn_x,
    syn_time = as.vector(cat_init$syn_time[, 1]),
    syn_status = as.vector(cat_init$syn_status[, 1]),
    hazard_constant = cat_init$hazard_constant,
    adj_risk_set = risk_and_failure_sets$risk_set - risk_and_failure_sets$failure_set,
    failure_set = risk_and_failure_sets$failure_set,
    alpha = alpha,
    gamma = gamma,
    kappa = get_cox_kappa(
      X = cat_init$syn_x,
      time = cat_init$syn_time[, 1],
      status = cat_init$syn_status[, 1],
      hazard_constant = cat_init$hazard_constant,
    )
  )

  # Stan model and sampling
  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    ...
  )

  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]
  names(coefs) <- colnames(cat_init$adj_x)

  cat_cox_bayes_joint <- list(
    function_name = "cat_cox_bayes_joint",
    formula = as.formula(paste0(
      "survival::Surv(",
      cat_init$time_col_name,
      ",",
      cat_init$status_col_name,
      ")",
      f_rhs
    )),
    cat_init = cat_init,
    stan_chains = stan_chains,
    stan_iter = stan_iter,
    stan_warmup = stan_warmup,
    c_0 = c_0,
    H_star_list = H_star_list,
    stan_model = stan_model,
    stan_sample_model = stan_sample_model,
    coefficients = coefs,
    increment_cumulative_baseline_hazard = rstan::summary(stan_sample_model)$summary[
      grep("h_j", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ],
    tau = rstan::summary(stan_sample_model)$summary[
      grep("tau", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_cox_bayes_joint) <- c(cat_cox_bayes_joint$class, "cat_cox_bayes_joint", "cat_cox_bayes", "bayes")


  return(cat_cox_bayes_joint)
}

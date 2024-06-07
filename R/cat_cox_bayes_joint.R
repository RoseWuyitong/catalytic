#' Title
#'
#' @param data_
#' @param ...
#' @param initial_values
#' @param hazard_constant
#' @param tol
#' @param max_iter
#' @param c_0 small value of c_0 corresponds to a more diffuse and less informative prior
#' @param alpha
#' @param gamma
#' @param chains
#' @param iter
#' @param warmup
#' @param H_star_list
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
cat_cox_bayes_joint <- function(formula,
                                data_,
                                ...,
                                alpha = 2,
                                gamma = 1,
                                chains = 4,
                                iter = 2000,
                                warmup = 1000,
                                initial_values = NULL,
                                hazard_constant = NULL,
                                tol = 1e-5,
                                max_iter = 25,
                                c_0 = 2,
                                H_star_list = NULL) {
  .not_true_stop(
    is.data_cox(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization_cox."
    )
  )

  if (is.null(hazard_constant)) hazard_constant <- sum(data_$obs_event) / mean(data_$obs_time) / data_$obs_size
  data_ <- get_adj_data(data_ = data_, formula = formula)

  # Get unique times from sorted observed event times and a value larger than the maximum observed time
  s <- unique(c(
    # Sort the times where events occurred
    sort(data_$obs_time[data_$obs_event == 1]),
    # Create a value larger than the maximum observed time
    2 * max(data_$obs_time) -
      max(data_$obs_time[-which(data_$obs_time == max(data_$obs_time))])
  ))

  # get risk set and failure set
  intervals <- get_interval_indicators(data_$obs_time, data_$obs_event, s)

  if (is.null(H_star_list)) {
    # Get the initial guess of H_0 from the intercept-only weibull regression model based on observation data
    surv_weibull_model <- do.call(
      survival::survreg,
      list(
        formula = survival::Surv(data_$obs_time, data_$obs_event) ~ 1,
        dist = "weibull"
      )
    )
    # Extracts the shape and scale parameters of the Weibull distribution
    weibull_shape <- 1 / surv_weibull_model$scale
    weibull_scale <- exp(coef(surv_weibull_model))
    # Calculates the cumulative hazard values for the intervals
    H_star_list <- (s / weibull_scale)^weibull_shape
  }

  # List of alpha_0 value
  alpha_0_list <- c_0 * H_star_list

  # List of alpha value from each h_j ~ Gamma(alpha_j, c_0), where alpha_j = alpha_0[j] - alpha_0[j-1]
  alpha_of_h_j_list <- diff(c(0, alpha_0_list))

  data <- list(
    J = length(s),
    alpha_of_h_j_list = alpha_of_h_j_list,
    c_0 = c_0,
    obs_size = data_$obs_size,
    col_num = ncol(data_$adj_x),
    obs_x = data_$adj_obs_x,
    syn_size = data_$syn_size,
    syn_x = data_$adj_syn_x,
    syn_time = as.vector(data_$syn_time),
    syn_event = as.vector(data_$syn_event),
    hazard_constant = data_$hazard_constant,
    adj_risk_set = intervals$risk_set - intervals$failure_set,
    failure_set = intervals$failure_set,
    alpha = alpha,
    gamma = gamma,
    kappa = get_cox_kappa(
      X = data_$syn_x,
      time = as.vector(data_$syn_time),
      event = as.vector(data_$event),
      hazard_constant = data_$hazard_constant,
    )
  )

  # Stan model and sampling
  sample_model <- rstan::sampling(
    object = stan_model_selection(
      method = "cox",
      joint_tau = TRUE
    ),
    data = data,
    chains = chains,
    iter = iter,
    warmup = warmup
  )

  cat_cox_obj <- list(
    function_name = "cat_cox_bayes_joint",
    formula = formula,
    data_ = data_,
    chains = chains,
    iter = iter,
    warmup = warmup,
    tol = tol,
    max_iter = max_iter,
    c_0 = c_0,
    H_star_list = H_star_list,
    sample_model = sample_model,
    coefficients = rstan::summary(sample_model)$summary[
      grep("coefs", rownames(rstan::summary(sample_model)$summary)),
      "mean"
    ],
    increment_cumulative_baseline_hazard = rstan::summary(sample_model)$summary[
      grep("h_j", rownames(rstan::summary(sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_cox_obj) <- c(cat_cox_obj$class, "bayes", "cat_glm_bayes")


  return(cat_cox_obj)
}

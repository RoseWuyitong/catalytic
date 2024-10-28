# ----- GLMs -----

#' Generate Sample Data for GLM
#'
#' This function generates sample data for a specified GLM family. It can generate
#' binomial or Gaussian distributed data based on the provided parameters.
#'
#' @param family Character. The family of the GLM. Options are `"binomial"` or `"gaussian"`.
#' @param n Integer. The number of samples to generate.
#' @param mean Numeric. The mean of the distribution (used for both binomial and Gaussian).
#' @param sd Numeric. The standard deviation of the distribution (used only for Gaussian).
#'
#' @return Numeric vector of generated sample data.
#'
#' @export
#'
#' @examples
#' # Generate sample data for binomial family
#' binomial_data <- get_glm_sample_data(family = "binomial", n = 10, mean = 0.5)
#' binomial_data
#'
#' # Generate sample data for Gaussian family
#' gaussian_data <- get_glm_sample_data(family = "gaussian", n = 10, mean = 0, sd = 1)
#' gaussian_data
get_glm_sample_data <- function(family_string,
                                n = 10,
                                mean = 0,
                                sd = 1) {
  switch(family_string,
    "binomial" = stats::rbinom(n, 1, mean),
    "gaussian" = stats::rnorm(n, mean, sd)
  )
}

get_hmc_mcmc_result <- function(
    neg_log_den_func,
    neg_log_den_grad_func,
    coefs_0,
    iter = 5,
    hmc_scale = 0.01,
    hmc_steps = 5) {
  mcmc_collect <- matrix(NA, iter + 1, length(coefs_0))
  mcmc_collect[1, ] <- coefs_0

  for (i in seq_len(iter)) {
    mcmc_collect[i + 1, ] <- hmc_neal_2010(
      neg_log_den_func = neg_log_den_func,
      neg_log_den_grad_func = neg_log_den_grad_func,
      leapfrog_stepsize = rbeta(1, 0.5, 0.5) * hmc_scale / hmc_steps, # add random value to leapfrog_stepsize to overcome poor performance
      leapfrog_step = hmc_steps,
      current_pos = mcmc_collect[i, ]
    )$position
  }


  return(mcmc_collect[iter + 1, ])
}


# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
#
# Radford M. Neal, 2010.
#
# This program appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# to appear in the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   neg_log_den_func          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   neg_log_den_grad_func     A function to evaluate the gradient of neg_log_den_func.
#
#   leapfrog_stepsize    The stepsize to use for the leapfrog steps.
#
#   leapfrog_step          The number of leapfrog steps to do to propose a new state.
#
#   current_pos  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_pos if the endpoint of the
# trajectory was rejected).
#
# This function was written for illustrative purposes.  More elaborate
# implementations of this basic HMC method and various variants of HMC
# are available from my web page, http://www.cs.utoronto.ca/~radford/
hmc_neal_2010 <- function(
    neg_log_den_func,
    neg_log_den_grad_func,
    leapfrog_stepsize,
    leapfrog_step,
    current_pos) {
  pos <- current_pos

  # Sample initial momentum from standard normal distribution
  momentum <- rnorm(length(pos), 0, 1)
  current_momentum <- momentum

  # Initial half-step for momentum
  momentum <- momentum - leapfrog_stepsize * neg_log_den_grad_func(pos) / 2

  # Full steps for position and momentum over the leapfrog trajectory (except at end of trajectory)
  for (i in seq_len(leapfrog_step - 1)) {
    # Full step for position
    pos <- pos + leapfrog_stepsize * momentum

    # Full step for momentum
    momentum <- momentum - leapfrog_stepsize * neg_log_den_grad_func(pos)
  }

  # Final position update
  pos <- pos + leapfrog_stepsize * momentum

  # Final half-step for momentum and reverse momentum to maintain proposal symmetry
  momentum <- -(momentum - leapfrog_stepsize * neg_log_den_grad_func(pos) / 2)

  # Compute the potential and kinetic energy at the start and end of the trajectory
  current_potential_energy <- neg_log_den_func(current_pos)
  current_kinetic_energy <- sum(current_momentum^2) / 2
  proposed_potential_energy <- neg_log_den_func(pos)
  proposed_kinetic_energy <- sum(momentum^2) / 2

  # Metropolis acceptance criterion
  if (runif(1) < exp(current_potential_energy - current_potential_energy +
    current_kinetic_energy - proposed_kinetic_energy)) {
    # Accept the state at end of trajectory, returning the position at the end of the trajectory
    return(list(
      position = pos,
      potential_energy = proposed_potential_energy,
      accepted = TRUE
    ))
  } else {
    # Reject the state at end of trajectory, returning the current position
    return(list(
      position = current_pos,
      potential_energy = current_potential_energy,
      accepted = FALSE
    ))
  }
}

#' Compute Gradient of Log Density for GLM Families
#'
#' This function calculates the gradient of the log density with respect to the coefs
#' for a given GLM family based on the provided predictors, response variable, and weights.
#'
#' @param family Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param Y Vector. The response variable.
#' @param coefs Numeric vector. The coefs for the GLM.
#' @param weights Numeric vector. The weights for the GLM. Default is 1.
#'
#' @return Numeric vector. The gradient of the log density with respect to the coefs.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' Y <- rbinom(5, 1, 0.5)
#' coefs <- rnorm(4)
#'
#' # Compute gradient of log density for a binomial family
#' grad_binomial <- get_glm_log_density_grad("binomial", X, Y, coefs)
#' grad_binomial
get_glm_log_density_grad <- function(family_string,
                                     X,
                                     Y,
                                     coefs,
                                     weights = 1) {
  # Calculate the gradient of log density
  gradient <- weights * (Y - get_glm_mean(family_string, X, coefs))

  return(c(t(gradient) %*% as.matrix(cbind(1, X))))
}

#' Compute Diagonal Approximate Covariance Matrix
#'
#' This function computes the diagonal elements of the approximate covariance matrix
#' for the coefs in a generalized linear model (GLM). The covariance is derived
#' from the second derivative (Hessian) of the log-likelihood function.
#'
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param model A fitted GLM model object. The object should contain the fitted values
#'        and prior weights necessary for computing the Hessian.
#'
#' @return Numeric vector. The diagonal elements of the approximate covariance matrix.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' y <- rbinom(5, 1, 0.5)
#' model <- glm(y ~ X, family = binomial)
#'
#' diag_cov <- get_glm_diag_approx_cov(X, model)
#' diag_cov
get_glm_diag_approx_cov <- function(X,
                                    model) {
  # Nested function to compute the second derivative (Hessian) matrix
  get_second_derivative <- function(X,
                                    fitted_values,
                                    prior_weights) {
    return((t(cbind(1, X))
    ) %*% diag(fitted_values * (1 - fitted_values) * prior_weights) %*% cbind(1, X))
  }

  # Calculate the diagonal elements of the covariance approximation
  hessian <- get_second_derivative(
    X = as.matrix(X),
    fitted_values = model$fitted.values,
    prior_weights = model$prior.weights
  )

  tryCatch(
    return(diag(solve(hessian))),
    error = function(.) {
      # Regularization for dealing with computationally singular issue
      return(diag(solve(hessian + diag(ncol(hessian)) * 1e-5)))
    }
  )
}



#' Generate Suggestions for Bayesian Joint GLM Parameter Estimation
#'
#' This function provides suggestions for improving the parameter estimation process
#' in Bayesian joint GLM modeling based on the diagnostic output from a Stan model.
#' It evaluates the results and suggests adjustments to improve model fit and convergence.
#'
#' @param family Character. The GLM family to use, e.g., `"binomial"`.
#' @param stan_iter Integer. The number of iterations used in Stan sampling.
#' @param alpha Numeric. The alpha parameter used in the prior distribution.
#' @param stan_sample_model Stan model object containing sampling results.
#' @param binomial_use_theta Logical. Whether to use theta in the binomial model.
#' @param binomial_joint_alpha Logical. Whether to use joint alpha in the binomial model.
#' @param binomial_tau_lower Numeric. The lower bound for tau in the binomial model.
#'
#' @return NULL. The function prints suggestions to the console based on the model diagnostics.
#'
#' @export
#'
#' @examples
#' # Assuming you have a Stan model and relevant parameters
#' # get_glm_bayes_joint_suggestion("binomial", 2000, 1, stan_model, TRUE, FALSE, 0.5)
print_glm_bayes_joint_suggestion <- function(
    alpha,
    stan_iter,
    stan_sample_model,
    binomial_use_theta,
    binomial_joint_alpha,
    binomial_tau_lower,
    gaussian_known_variance) {
  # Helper function to calculate the error between estimated and theoretical percentiles
  get_est_params <- function(params, percentiles, mean_, distribution) {
    alpha <- exp(params[1]) # Convert log-alpha to alpha
    beta <- exp(params[2]) # Convert log-beta to beta

    # Calculate percentiles based on the distribution
    calc_percentiles <- switch(distribution,
      "invgamma" = invgamma::qinvgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta),
      "gamma" = stats::qgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta)
    )

    # Compute error between calculated percentiles and those from the Stan model
    error <- sum((percentiles - calc_percentiles)^2) + (mean_ - (beta / (alpha - 1)))^2

    return(error)
  }

  if (binomial_use_theta) {
    # Retrieve theta summary from Stan model
    theta_summary <- rstan::summary(stan_sample_model)$summary["theta", ]

    # Optimize parameters for the inverse gamma distribution
    est_params <- stats::optim(
      par = c(1, 1),
      fn = get_est_params,
      percentiles = theta_summary[c("25%", "50%", "75%")],
      mean_ = theta_summary["mean"],
      distribution = "invgamma"
    )

    # Check if estimated parameters fall outside the 97.5% quantile range
    if (invgamma::qinvgamma(c(0.975),
      shape = exp(est_params$par[1]),
      scale = exp(est_params$par[2])
    ) < theta_summary["97.5%"]) {
      # Provide suggestions based on model diagnostics
      if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
        cat("\nSUGGESTION: Please consider using more iterations for Stan (e.g.,", stan_iter + 1000, ")\n")
      } else if (alpha < 5 & !binomial_joint_alpha) {
        cat("\nSUGGESTION: Please consider using a larger `alpha` (e.g.,", min(alpha + 0.1, 5), ")\n")
      } else if (alpha >= 5 & !binomial_joint_alpha) {
        cat("\nSUGGESTION: Please consider setting `binomial_joint_alpha = TRUE`.\n")
      }
    }
  } else {
    # Retrieve tau summary from Stan model
    tau_summary <- rstan::summary(stan_sample_model)$summary["tau", ]

    # Optimize parameters for the gamma distribution
    est_params <- stats::optim(
      par = c(1, 1),
      fn = get_est_params,
      percentiles = tau_summary[c("25%", "50%", "75%")],
      mean_ = tau_summary["mean"],
      distribution = "gamma"
    )

    # Check if estimated parameters fall outside the 2.5% quantile range
    if (stats::qgamma(c(0.025),
      shape = exp(est_params$par[1]),
      scale = exp(est_params$par[2])
    ) > tau_summary["2.5%"]) {
      # Provide suggestions based on model diagnostics
      if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
        cat("\nSUGGESTION: Please consider using more iterations for Stan (e.g.,", stan_iter + 1000, ")\n")
      } else if (binomial_tau_lower < 1) {
        cat("\nSUGGESTION: Please consider using a larger `binomial_tau_lower` (e.g.,", min(binomial_tau_lower + 0.05, 1), ")\n")
      } else {
        cat("\nSUGGESTION: Please consider setting `binomial_use_theta = TRUE`.\n")
      }
    }
  }
}


get_glm_custom_var <- function(
    cat_init,
    tau,
    coefs = NULL,
    formula = NULL) {
  if (is.null(coefs)) {
    coefs <- coef(
      do.call(
        stats::glm,
        list(
          formula = formula,
          family = cat_init$family,
          data = cat_init$data,
          weights = c(
            rep(1, cat_init$obs_size),
            rep(tau / cat_init$syn_size, cat_init$syn_size)
          )
        )
      )
    )
  }

  return((get_discrepancy(
    discrepancy_method = "square_error",
    family_string = cat_init$family,
    Y = cat_init$obs_y,
    X = if ("adj_obs_x" %in% cat_init) cat_init$adj_obs_x else cat_init$obs_x,
    coefs = coefs
  ) + (tau / cat_init$syn_size) * get_discrepancy(
    discrepancy_method = "square_error",
    family_string = cat_init$family,
    Y = cat_init$syn_y,
    X = if ("adj_syn_x" %in% cat_init) cat_init$adj_syn_x else cat_init$syn_x,
    coefs = coefs
  )
  ) / ((cat_init$obs_size + tau - ncol(cat_init$adj_x) + 2)))
}

#' Compute Lambda Based on Discrepancy Method
#'
#' This function calculates a lambda value based on the selected discrepancy method
#' for a generalized linear model (GLM). The discrepancy method determines the
#' type of error or deviance used in the calculation.
#'
#' @param discrepancy_method Character. A string specifying the type of discrepancy method
#'        to use. Options are `"square_error"`, `"classification_error"`, or `"logistic_deviance"`.
#'        Default is `"square_error"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param coefs Numeric vector. The coefficients for the GLM.
#'
#' @return Numeric. The computed lambda value based on the selected discrepancy method.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' coefs <- rnorm(4)
#'
#' lambda <- get_glm_lambda(discrepancy_method = "square_error", X, coefs)
#' lambda
get_glm_lambda <- function(
    discrepancy_method = c(
      "square_error",
      "classification_error",
      "logistic_deviance"
    ),
    X,
    coefs) {
  # Validate and match the discrepancy method
  discrepancy_method <- match.arg(discrepancy_method)

  # Compute lambda based on the selected discrepancy method
  return(switch(discrepancy_method,
    "square_error" = 2 * get_glm_mean(family_string = "gaussian", X = X, coefs = coefs),
    "classification_error" = 2 * get_glm_mean(family_string = "binomial", X = X, coefs = coefs),
    "logistic_deviance" = get_linear_predictor(X = X, coefs = coefs)
  ))
}

#' Compute Log Density Based on GLM Family
#'
#' This function calculates the log density of the response variable given a generalized linear model (GLM)
#' based on the specified family. The log density is computed differently for binomial and gaussian families.
#'
#' @param family Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param Y Vector or data frame. The response variable for the GLM. If a data frame, it is converted to a numeric vector.
#' @param coefs Numeric vector. The coefficients for the GLM.
#' @param weights Numeric vector. Weights for the observations. Default is `1` (no weighting).
#'
#' @return Numeric. The computed log density of the response variable based on the specified family.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' Y <- rbinom(5, 1, 0.5)
#' coefs <- rnorm(4)
#'
#' # Compute log density for a binomial family
#' log_density_binomial <- get_glm_log_density("binomial", X, Y, coefs)
#' log_density_binomial
get_glm_log_density <- function(family_string,
                                X,
                                Y,
                                coefs,
                                weights = 1) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X = X, coefs = coefs)

  return((weights * switch(family_string,
    "binomial" = sum((Y * lp - log(1 + exp(lp)))),
    "gaussian" = -0.5 * sum((log(2 * pi * max(mean((Y - lp)^2), 1e-5)) + 1)) # TODO: confirm w/ prof
  )))
}


#' Compute Mean Based on GLM Family
#'
#' This function calculates the mean of the response variable for a generalized linear model (GLM)
#' based on the specified family. The calculation depends on whether the family is binomial or gaussian.
#'
#' @param family Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param coefs Numeric vector. The coefficients for the GLM.
#'
#' @return Numeric vector. The computed mean of the response variable based on the specified family.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' coefs <- rnorm(4)
#'
#' # Compute mean for a binomial family
#' mean_binomial <- get_glm_mean("binomial", X, coefs)
#' mean_binomial
#'
#' # Compute mean for a gaussian family
#' mean_gaussian <- get_glm_mean("gaussian", X, coefs)
#' mean_gaussian
get_glm_mean <- function(family_string, X, coefs) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X = X, coefs = coefs)

  # Return the mean based on the specified family
  return(switch(family_string,
    "binomial" = 1 / (1 + exp(-lp)),
    "gaussian" = lp
  ))
}

validate_glm_input <- function(formula,
                               cat_init,
                               tau = NULL,
                               tau_seq = NULL,
                               tau_0 = NULL,
                               parametric_bootstrap_iteration_times = NULL,
                               cross_validation_fold_num = NULL,
                               risk_estimate_method = NULL,
                               discrepancy_method = NULL,
                               binomial_tau_lower = NULL,
                               tau_alpha = NULL,
                               tau_gamma = NULL,
                               gibbs_iter = NULL,
                               gibbs_warmup = NULL,
                               coefs_iter = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau, incl_0 = TRUE)
  validate_positive("tau_0", tau_0, incl_0 = TRUE)
  validate_positive("parametric_bootstrap_iteration_times", parametric_bootstrap_iteration_times)
  validate_positive("cross_validation_fold_num", cross_validation_fold_num)
  validate_positive("binomial_tau_lower", binomial_tau_lower)
  validate_positive("tau_alpha", tau_alpha)
  validate_positive("tau_gamma", tau_gamma)
  validate_positive("gibbs_iter", gibbs_iter)
  validate_positive("gibbs_warmup", gibbs_warmup)
  validate_positive("coefs_iter", coefs_iter)
  validate_positive("tau_seq", tau_seq, incl_0 = TRUE, is_vector = TRUE)



  # Check if cat_init is generated from cat_glm_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_glm_initialization")) {
    stop(
      paste0(
        "The provided `cat_init` must be generated from `cat_glm_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if formula has the same response name as that from cat_init
  f_lhs <- get_formula_lhs(formula)
  if (f_lhs != "NULL" & f_lhs != cat_init$y_col_name) {
    stop(
      paste0(
        "The response name from provided `formula` <",
        get_formula_string(formula),
        "> should be same as the response name used in the provided `cat_init` object."
      ),
      call. = FALSE
    )
  }

  # Check for gibbs
  if (!is.null(gibbs_iter) & !is.null(gibbs_warmup)) {
    # Check if user apply gibbs to linear regression
    if (cat_init$family != "binomial") {
      stop(
        "\nPlease note that only logistic regression can use `cat_glm_bayes_joint_gibbs`.",
        call. = FALSE
      )
    }

    # Check if warmup is larger than iter
    if (gibbs_warmup >= gibbs_iter) {
      stop(
        paste0(
          "The provided `warmup` should be less than `iter`. <",
          "iter = ",
          gibbs_iter,
          "; warmup =",
          gibbs_warmup,
          "> is invalid.",
        ),
        call. = FALSE
      )
    }

    # Check if coefs_iter too large in gibbs
    if (coefs_iter > 10) {
      response <- menu(
        c("YES", "NO"),
        title = paste0(
          "\nPlease note that large `coefs_iter` can slow down gibb's sampling process.\n",
          "Do you want to continue?\n"
        )
      )

      if (response == 2) {
        stop("\nProcess terminated by user.", call. = FALSE)
      }
    }
  }

  # Specific check for cat_glm_tune
  if (!is.null(risk_estimate_method) & !is.null(discrepancy_method)) {
    # Initialize an empty vector for validation errors
    validation_errors <- c()

    # Append messages based on conditions
    if (cat_init$family == "gaussian" & risk_estimate_method == "steinian_estimate") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "\nOnly logistic regression can use `risk_estimate_method = steinian_estimate`.\n",
          "SUGGESTION: Use `mallowian_estimate` for linear regression."
        )
      )
    }

    if (cat_init$family == "binomial" & risk_estimate_method == "mallowian_estimate") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "\nOnly linear regression can use `risk_estimate_method = mallowian_estimate`.\n",
          "SUGGESTION: Use `steinian_estimate` for logistic regression."
        )
      )
    }

    if (cat_init$family == "gaussian" & discrepancy_method != "square_error") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "\nOnly logistic regression can use `discrepancy_method = classification_error` or `logistic_deviance`.\n",
          "SUGGESTION: Use `square_error` for linear regression."
        )
      )
    }

    if (cat_init$family == "binomial" & discrepancy_method == "square_error") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "\nOnly linear regression can use `discrepancy_method = square_error`.\n",
          "SUGGESTION: Use `logistic_deviance` for logistic regression."
        )
      )
    }

    # Stop the process if there are any validation errors
    if (length(validation_errors) > 0) {
      stop(paste(validation_errors, collapse = "\n"))
    }

    # Warn User that
    if (cat_init$family == "binomial" &
      risk_estimate_method == "steinian_estimate" &
      discrepancy_method != "logistic_deviance") {
      warning(
        paste0(
          "\nOnly `discrepancy_method = logistic_deviance` can be applied to `risk_estimate_method = steinian_estimate`.\n",
          "Will continue use `discrepancy_method = logistic_deviance` instead of <",
          discrepancy_method,
          ">."
        )
      )
    }


    # Set data size thresholds based on the risk estimation method
    size_threshold <- switch(risk_estimate_method,
      "parametric_bootstrap" = 1000,
      "steinian_estimate" = 800,
      Inf
    )

    # Check if the observational dataset size exceeds the threshold
    if (cat_init$obs_size > size_threshold) {
      response <- menu(
        choices = c("YES", "NO"),
        title = paste0(
          "The large size of the dataset may lead to longer estimation times.\n",
          "Do you want to continue?"
        )
      )


      # Terminate process if the user chooses "NO"
      if (response == 2) {
        stop("\nProcess terminated by user.", call. = FALSE)
      }
    }

    # Issue warning for larger dataset sizes
    if (cat_init$obs_size > (size_threshold / 2)) {
      warning(
        paste0(
          "The large size of the dataset may lead to longer estimation times.\n",
          "SUGGESTION: Consider using `risk_estimate_method = cross_validation` ",
          "for quicker estimation."
        ),
        call. = FALSE
      )
    }

    # Check if dataset size is small for cross-validation
    if (risk_estimate_method == "cross_validation" & cat_init$obs_size < (5 * cross_validation_fold_num)) {
      warning(
        paste0(
          "The size of the observation dataset is relatively small for cross-validation to yield optimal results.\n",
          "SUGGESTION: Please consider utilizing alternative methods (e.g., risk_estimate_method = parametric_bootstrap) ",
          "or reduce `cross_validation_fold_num` for better estimation."
        ),
        call. = FALSE
      )
    }
  }
}


#' Retrieve GLM Family Name or Name with Link Function
#'
#' This function retrieves the name of a GLM family or, optionally, the family name with the associated link function.
#'
#' @param family Character or function. The name of the GLM family (as a string) or a function that returns a GLM family object.
#' @param with_link Logical. If TRUE, returns the family name along with the link function in the format "family [link]". If FALSE, only the family name is returned. Default is FALSE.
#'
#' @return A character string. The name of the GLM family, or the name with the link function if `with_link` is TRUE.
#'
#' @export
#'
#' @examples
#' # Retrieve the name of the Gaussian family
#' family_name <- get_glm_family_string("gaussian")
#' family_name
#'
#' # Retrieve the name of the Poisson family with its link function
#' family_name_with_link <- get_glm_family_string("poisson", with_link = TRUE)
#' family_name_with_link
get_glm_family_string <- function(family, with_link = FALSE) {
  # Validate and retrieve the family object
  validate_glm_family <- function(f) {
    if (is.character(f)) {
      tryCatch(
        {
          f <- get(f, mode = "function", envir = parent.frame(2))
        },
        error = function(e) {
          stop("'family' must be a valid GLM family.",
            call. = FALSE
          )
        }
      )
    }
    if (is.function(f)) {
      f <- f()
    }
    if (!is(f, "family")) {
      stop("'family' must be a valid GLM family.",
        call. = FALSE
      )
    }
    return(f)
  }

  # Retrieve and validate the family object
  family <- validate_glm_family(family)

  # Return the family name or family name with link function
  if (with_link) {
    return(paste0(family$family, " [", family$link, "]"))
  } else {
    return(family$family)
  }
}

validate_glm_initialization_input <- function(formula,
                                              family,
                                              data,
                                              syn_size,
                                              custom_variance,
                                              x_degree) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)
  validate_positive("custom_variance", custom_variance)
  validate_positive("x_degree", x_degree, is_vector = TRUE)

  # Check if formula is for GLMs
  if (grepl("\\(.*\\|.*\\)", deparse(formula)) || grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for Generalized Linear Models (GLMs)."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if the formula is overly complex
  if (length(attr(terms(formula), "term.labels")) >= (ncol(data) - 1)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }

  # Check if family is either binomial or gaussian
  family_string <- get_glm_family_string(family)
  if (!(family_string %in% c("binomial", "gaussian"))) {
    stop(
      paste0(
        "The provided `family` must be either `binomial` or `gaussian`. <",
        family,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if x_degree has same size to the data covariance
  if (!is.null(x_degree) & (length(x_degree) != (ncol(data) - 1))) {
    stop(
      paste0(
        "The provided `x_degree` should be a vector of positive values ",
        "with a length less than the number of covariates. <",
        paste(x_degree, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  if (!is.null(syn_size) & syn_size < (4 * ncol(data))) {
    warning(
      "The provided `syn_size` should is too small, recommand to choose a larger value.",
      call. = FALSE
    )
  }

  if (family_string != "gaussian" & (!is.null(custom_variance))) {
    warning(
      "The `custom_variance` is not applicable when the family is not `gaussian`.",
      call. = FALSE
    )
  }
}

# ----- COX -----

#' Calculate Risk and Failure Sets for Cox Proportional Hazards Model
#'
#' This function calculates the risk and failure sets for subjects in a Cox proportional hazards model based on their time-to-event data, status, and an indicator vector.
#'
#' @param time_vector A numeric vector of time-to-event data for each subject.
#' @param status_vector A numeric vector indicating event occurrence (1 = event, 0 = censored).
#' @param indicator_vector A numeric vector representing the indicator times used to define risk and failure sets.
#'
#' @return A list containing two elements:
#' \item{risk_set}{A matrix indicating which subjects are at risk at each time point.}
#' \item{failure_set}{A matrix indicating which subjects experienced an event at each time point.}
#'
#' @export
#'
#' @examples
#' time_vector <- rexp(10)
#' status_vector <- rbinom(10, 1, 0.5)
#' indicator_vector <- sort(unique(c(time_vector, time_vector + rexp(10))))
#'
#' # Calculate risk and failure sets
#' sets <- get_cox_risk_and_failure_sets(time_vector, status_vector, indicator_vector)
#' sets$risk_set
#' sets$failure_set
get_cox_risk_and_failure_sets <- function(time_vector,
                                          status_vector,
                                          indicator_vector) {
  time_vector <- as.vector(time_vector)
  status_vector <- as.vector(status_vector)
  indicator_vector <- as.vector(indicator_vector)

  # Function to find the index where an event occurs or reaches censoring time
  find_occur_idx <- function(time) {
    min(which(indicator_vector >= time))
  }


  # Initialize risk and failure sets as zero matrices
  risk_set <- failure_set <- matrix(
    0,
    length(time_vector),
    length(indicator_vector)
  )

  # Process subjects who experienced an event (status = 1)
  for (i in which(status_vector == 1)) {
    # Mark the exact time of the event in the failure set
    failure_set[i, find_occur_idx(time_vector[i])] <- 1
    # Mark all times up to the event time in the risk set
    risk_set[i, 1:find_occur_idx(time_vector[i])] <- 1
  }

  # Mark all times up to the censoring time in the risk set
  risk_set[which(status_vector == 0), 1:find_occur_idx(time_vector[i])] <- 1

  # Return the risk and failure sets as a list
  return(list(
    risk_set = risk_set,
    failure_set = failure_set
  ))
}


#' Estimate the kappa value for the synthetic Cox proportional hazards model
#'
#' This function iteratively estimates the kappa value for the synthetic Cox proportional hazards model using a vectorized approach for efficiency.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param status A vector indicating event occurrence (1 = event, 0 = censored).
#' @param hazard_constant A scalar representing the hazard constant. Defaults to NULL, in which case it's calculated internally.
#' @param tol Tolerance for convergence. Defaults to 1e-5.
#' @param max_iter Maximum number of iterations for the estimation process. Defaults to 25.
#' @param entry_points Optional vector indicating entry times for the observations. Defaults to NULL.
#' @param initial_values Optional vector of initial coefficient values. Defaults to NULL.
#'
#' @return A numeric value representing the estimated kappa for the synthetic Cox model.
#'
#' @export
#'
#' @examples
#' # Randomly generated data for demonstration
#' X <- matrix(rnorm(100), ncol = 5)
#' time <- rexp(20)
#' status <- rbinom(20, 1, 0.5)
#' hazard_constant <- 0.02
#'
#' # Estimate kappa
#' kappa <- get_cox_kappa(X, time, status, hazard_constant)
#' kappa
get_cox_kappa <- function(X,
                          time,
                          status,
                          hazard_constant,
                          entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize coefficients
  coefs <- rep(0, ncol(X))

  diff_new_old <- TRUE
  iter <- 0

  while ((iter < 25) & diff_new_old) {
    coefs_old <- coefs

    # Calculate linear predictor
    exp_lp <- exp(c(get_linear_predictor(X, coefs_old)))

    # Calculate gradient and Hessian
    gd_old <- get_cox_syn_gradient(X, time, coefs_old, hazard_constant)
    hess_old <- get_cox_syn_hessian(X, time, coefs_old, hazard_constant)

    # Update coefficients using QR solve
    coefs <- coefs_old - get_cox_qr_solve(hess_old, gd_old)

    iter <- iter + 1

    # Check for convergence
    diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > -5 * norm(as.numeric(coefs_old), "2")
  }

  # Calculate the final linear predictor
  lp_new <- c(get_linear_predictor(X, coefs))

  return(sum(lp_new - exp(lp_new) * hazard_constant * time) / nrow(X))
}

#' Compute the Partial Likelihood for the Cox Proportional Hazards Model
#'
#' This function calculates the partial likelihood for the Cox proportional hazards model. The partial likelihood is computed for the censored observations in the dataset.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param status A vector indicating the event status (1 for event occurred, 0 for censored).
#' @param coefs A vector of regression coefficients.
#' @param entry_points A vector of entry points (optional). Defaults to NULL, in which case a vector of zeros is used.
#'
#' @return A numeric scalar representing the partial likelihood of the Cox model.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(20), ncol = 4)
#' time <- runif(5, 0, 10)
#' status <- sample(0:1, 5, replace = TRUE)
#' coefs <- rnorm(4)
#' entry_points <- rep(0, 5)
#'
#' # Compute partial likelihood
#' pl <- get_cox_partial_likelihood(X, time, status, coefs, entry_points)
#' pl
get_cox_partial_likelihood <- function(X,
                                       time,
                                       status,
                                       coefs,
                                       entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  pl <- 0

  # Calculate partial likelihood for each censored observation
  for (i in which(status == 1)) {
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    # Compute the linear predictors for the risk set
    exp_risk_lp <- exp(
      c(get_linear_predictor(
        X[risk_set_idx, , drop = FALSE],
        coefs
      ))
    )

    # Update partial likelihood
    pl <- pl + sum(X[i, , drop = FALSE] * coefs) - log(sum(exp_risk_lp))
  }

  return(pl)
}


#' Solve Linear System using QR Decomposition
#'
#' @param coef_matrix A matrix of coefficients representing the system of linear equations.
#' @param constant_vector A numeric vector representing the constants in the linear system.
#'
#' @return A numeric vector representing the solution to the linear system. NA values in the solution
#' are replaced with a small value (0.0001). If an error occurs during solving, a vector of default values
#' (0.0001) is returned.
#'
#' @export
#'
#' @details This function solves the linear system defined by \code{coef_matrix} and \code{constant_vector}
#' using QR decomposition. Any NA values in the resulting solution vector are replaced with 0.0001. If there
#' is an error during the solution process, a vector of default values (0.0001) is returned instead.
#'
#' @examples
#' # Define a matrix of coefficients and a constant vector
#' coef_matrix <- matrix(c(1, 2, 3, 4), ncol = 2)
#' constant_vector <- c(5, 6)
#'
#' # Solve the linear system
#' solution <- get_cox_qr_solve(coef_matrix, constant_vector)
#' solution
#'
#' # Example with a matrix that might produce NA values
#' coef_matrix_na <- matrix(c(1, 1, 1, 1), ncol = 2)
#' constant_vector_na <- c(1, 1)
#'
#' # Solve the linear system and check the result
#' solution_na <- get_cox_qr_solve(coef_matrix_na, constant_vector_na)
#' solution_na
get_cox_qr_solve <- function(hessian_matrix,
                             gradient_vector) {
  # Calculate solution using QR decomposition with error handling
  qr_solution <- tryCatch(
    {
      # Attempt to solve the linear system
      solve(hessian_matrix, gradient_vector)
    },
    error = function(e) {
      # In case of error, return a vector of default values
      rep(0.0001, length(gradient_vector))
    }
  )

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
}


#' Compute the Gradient for Cox Proportional Hazards Model
#'
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param status A numeric vector of event indicators (1 for event, 0 for censored).
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param entry_points (Optional) A numeric vector of entry times for the subjects. Defaults to 0.
#'
#' @return A numeric vector representing the gradient of the Cox proportional hazards model.
#'
#' @export
#'
#' @details This function computes the gradient for the Cox proportional hazards model. The gradient
#' is calculated by considering the contributions of each observation to the gradient based on the
#' risk set at each event time.
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' time <- rexp(5, rate = 0.1)
#' status <- rbinom(5, 1, 0.5)
#' coefs <- rnorm(3)
#'
#' # Compute the gradient
#' gradient <- get_cox_gradient(X, time, status, coefs)
#' gradient
get_cox_gradient <- function(X,
                             time,
                             status,
                             coefs,
                             entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize gradient matrix
  gradient <- matrix(NA, nrow = ncol(X), ncol = nrow(X))

  # Loop over each observation to compute the gradient
  for (i in seq_len(nrow(X))) {
    # Get the indices of the risk set for the current observation
    risk_set_idx <- get_cox_risk_set_idx(
      time[i],
      entry_points,
      time,
      status
    )

    # Skip if the risk set is empty
    if (length(risk_set_idx) == 0) next

    # Extract risk set covariates and compute exponential risk
    risk_set <- as.matrix(X[risk_set_idx, , drop = FALSE])
    exp_risk_set <- exp(risk_set %*% coefs)

    # Compute the gradient for the current observation
    gradient[, i] <- c(X[i, , drop = FALSE]) -
      c(t(risk_set) %*% exp_risk_set / sum(exp_risk_set))
  }

  # Return the gradient vector by summing over all observations weighted by the status vector
  return(c(gradient %*% status))
}

#' Compute the gradient of the synthetic Cox proportional hazards model
#'
#' This function calculates the gradient of the synthetic Cox proportional hazards model using a vectorized approach.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param coefs A vector of regression coefficients.
#' @param hazard_constant A scalar representing the hazard constant.
#'
#' @return A numeric vector representing the gradient of the synthetic Cox model.
#'
#' @export
#'
#' @examples
#' X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' time <- c(1, 2, 3)
#' coefs <- c(0.5, -0.3)
#' hazard_constant <- 0.7
#'
#' # Compute gradient
#' gradient <- get_cox_syn_gradient(X, time, coefs, hazard_constant)
#' gradient
get_cox_syn_gradient <- function(X,
                                 time,
                                 coefs,
                                 hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate linear predictors
  syn_lp <- c(get_linear_predictor(X, coefs))

  # Calculate the gradient in a vectorized manner
  gradient <- c(colSums((1 - hazard_constant * time * exp(syn_lp)) * X))

  return(gradient)
}


#' Identify the risk set indices for Cox proportional hazards model
#'
#' This function returns the indices of the risk set for a given time of interest in the Cox proportional hazards model.
#'
#' @param time_of_interest A numeric value representing the time at which the risk set is calculated.
#' @param entry_vector A numeric vector representing the entry times of subjects.
#' @param time_vector A numeric vector representing the time-to-event or censoring times of subjects.
#' @param status_vector A numeric vector indicating event occurrence (1) or censoring (0) for each subject.
#'
#' @return A vector of indices representing the subjects at risk at the specified time of interest.
#'
#' @export
#'
#' @examples
#' entry_vector <- c(0, 0, 0, 0, 0)
#' time_vector <- c(3, 5, 8, 9, 10)
#' status_vector <- c(1, 1, 0, 1, 0)
#' time_of_interest <- 5
#'
#' # Identify risk set indices
#' risk_set_idx <- get_cox_risk_set_idx(time_of_interest, entry_vector, time_vector, status_vector)
#' risk_set_idx
get_cox_risk_set_idx <- function(time_of_interest,
                                 entry_vector,
                                 time_vector,
                                 status_vector) {
  time_vector <- as.vector(time_vector)
  entry_vector <- as.vector(entry_vector)
  status_vector <- as.vector(status_vector)

  # Find indices where subjects are at risk at the given time of interest
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & status_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}


#' Compute the Hessian Matrix for Cox Proportional Hazards Model
#'
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param status A numeric vector of event indicators (1 for event, 0 for censored).
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param entry_points (Optional) A numeric vector of entry times for the subjects. Defaults to 0.
#'
#' @return A matrix representing the negative Hessian of the Cox model.
#'
#' @export
#'
#' @details This function computes the Hessian matrix of the Cox proportional hazards model,
#' which is used for estimating the covariance matrix of the coefficients. The Hessian is calculated
#' by summing contributions from each event time in the risk set.
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' time <- rexp(5, rate = 0.1)
#' status <- rbinom(5, 1, 0.5)
#' coefs <- rnorm(3)
#'
#' # Compute the Hessian matrix
#' hessian_matrix <- get_cox_hessian(X, time, status, coefs)
#' hessian_matrix
get_cox_hessian <- function(X,
                            time,
                            status,
                            coefs,
                            entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize the Hessian matrix
  hessian <- matrix(0, nrow = ncol(X), ncol = ncol(X))

  # Calculate the exponential of the linear predictor
  exp_risk_lp <- exp(c(get_linear_predictor(X, coefs)))

  # Compute the Hessian matrix contributions for each censored event
  for (i in which(status == 1)) {
    # Determine the risk set for the current time of interest
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    # Skip if the risk set is empty
    if (length(risk_set_idx) == 0) next

    # Extract the risk set and compute relevant quantities
    risk_set <- as.matrix(X[risk_set_idx, , drop = FALSE])
    exp_risk_lp_risk_set <- exp_risk_lp[risk_set_idx]
    sum_exp_risk_lp <- sum(exp_risk_lp_risk_set)

    # Compute the Hessian matrix contribution
    risk_exp_risk_lp <- t(risk_set) %*% exp_risk_lp_risk_set
    risk_sqrt_exp_risk_lp <- t(sweep(risk_set, 1, sqrt(exp_risk_lp_risk_set), "*"))

    hessian <- hessian +
      tcrossprod(risk_sqrt_exp_risk_lp) / sum_exp_risk_lp -
      tcrossprod(risk_exp_risk_lp) / sum_exp_risk_lp^2
  }

  return(-hessian)
}

#' Compute the Synthetic Hessian Matrix for Cox Proportional Hazards Model
#'
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param hazard_constant A numeric value representing the hazard constant.
#'
#' @return A matrix representing the synthetic Hessian of the Cox model.
#'
#' @export
#'
#' @details This function computes the synthetic Hessian matrix for the Cox proportional hazards model.
#' The Hessian is calculated by summing the contributions from each individual observation, scaled by
#' the hazard constant and the time of the event.
#'
#' @examples
#' X <- matrix(rnorm(15), ncol = 3)
#' time <- rexp(5, rate = 0.1)
#' coefs <- rnorm(3)
#' hazard_constant <- 1
#'
#' # Compute the synthetic Hessian matrix
#' hessian_matrix <- get_cox_syn_hessian(X, time, coefs, hazard_constant)
#' hessian_matrix
get_cox_syn_hessian <- function(X,
                                time,
                                coefs,
                                hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate the linear predictors
  syn_lp <- c(get_linear_predictor(X, coefs))

  # Compute the Hessian using a vectorized approach
  hessian <- hazard_constant *
    t(X) %*% (sweep(X, 1, time * exp(syn_lp), "*"))

  return(-hessian)
}

validate_cox_input <- function(formula,
                               cat_init,
                               tau = NULL,
                               tau_seq = NULL,
                               init_coefficients = NULL,
                               tol = NULL,
                               max_iter = NULL,
                               cross_validation_fold_num = NULL,
                               hazard_beta = NULL,
                               tau_alpha = NULL,
                               tau_gamma = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau, incl_0 = TRUE)
  validate_positive("tol", tol)
  validate_positive("max_iter", max_iter)
  validate_positive("cross_validation_fold_num", cross_validation_fold_num)
  validate_positive("hazard_beta", hazard_beta)
  validate_positive("tau_alpha", tau_alpha)
  validate_positive("tau_gamma", tau_gamma)
  validate_positive("tau_seq", tau_seq, incl_0 = TRUE, is_vector = TRUE)

  # Check if cat_init is generated from cat_cox_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_cox_initialization")) {
    stop(
      paste0(
        "The provided `cat_init` must be generated from `cat_cox_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if formula has the same time and status column name as that from cat_init
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
        ",\\s*([^\\)]+)\\)",
        f_lhs
      )
    )[[1]][2]
  )

  if (f_lhs != "NULL" & (time_col_name != cat_init$time_col_name | status_col_name != cat_init$status_col_name)) {
    stop(
      paste0(
        "The time column or/and status column name from provided `formula` <",
        get_formula_string(formula),
        "> should be same as the response name used in the provided `cat_init` object."
      ),
      call. = FALSE
    )
  }


  # Check if the length of input coefs is match to the data.
  if (!is.null(init_coefficients) && length(init_coefficients) != ncol(cat_init$adj_x)) {
    stop(
      paste0(
        "The provided `init_coefficients` should be same length of number of columns from data."
      ),
      call. = FALSE
    )
  }

  # Check if dataset size is small for cross-validation
  if (!is.null(cross_validation_fold_num) && (cat_init$obs_size < (5 * cross_validation_fold_num))) {
    warning(
      paste0(
        "The size of the observational dataset is relatively small for cross-validation to yield results.\n",
        "SUGGESTION: Please consider reducing `cross_validation_fold_num` for better estimation."
      ),
    )
  }
}


validate_cox_initialization_input <- function(formula,
                                              data,
                                              syn_size,
                                              hazard_constant,
                                              entry_points,
                                              x_degree) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)
  validate_positive("hazard_constant", hazard_constant)
  validate_positive("entry_points", entry_points, is_vector = TRUE)

  # Check if formula is for COX
  if (!grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for a Cox Regression Models (COX)."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "> is invalid."
      ),
      call. = FALSE
    )
  }


  # Check if the formula is overly complex
  if (length(attr(terms(formula), "term.labels")) >= (ncol(data) - 1)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }
  # Check if x_degree is positive and same size to the data covariance
  if (!is.null(x_degree) && (length(x_degree) != (ncol(data) - 2) || any(x_degree < 0))) {
    stop(
      paste0(
        "The provided `x_degree` should be a vector of non-negtive values ",
        "with a length same as the number of covariates. <",
        paste(x_degree, collapse = ","),
        "> is invalid."
      )
    )
  }

  # Check if entry_points is positive and same size to the data covariance
  if (!is.null(entry_points) & (length(entry_points) != nrow(data))) {
    stop(
      paste0(
        "The provided `entry_points` should be a vector of values ",
        "with a length same as the number of covariates. <",
        paste(entry_points, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }
}
# ----- LMM -----

#' Update the log marginal likelihood for sigma^2 in a hierarchical linear mixed model with catalytic prior
#'
#' This function calculates the log marginal likelihood for sigma^2 in a hierarchical linear mixed model (LMM)
#' with a catalytic prior. The marginal likelihood is updated based on provided parameters and data sequences.
#'
#' @param residual_variance A scalar representing the variance of the residual errors.
#' @param random_effect_variance A scalar representing the variance of the random effects.
#' @param d_seq A numeric vector representing the diagonal elements of the design matrix.
#' @param d_seq_star A numeric vector representing the diagonal elements of the design matrix for the prior.
#' @param q A numeric vector representing the mean of the posterior distribution for the random effects.
#' @param q_star A numeric vector representing the mean of the prior distribution for the random effects.
#' @param tau A scalar representing the precision parameter of the catalytic prior.
#'
#' @return A scalar representing the updated log marginal likelihood for sigma^2.
#'
#' @export
#'
#' @examples
#' residual_variance <- 1
#' random_effect_variance <- 0.5
#' d_seq <- c(2, 3)
#' d_seq_star <- c(1, 2)
#' q <- c(0.5, 1)
#' q_star <- c(0.3, 0.7)
#' tau <- 0.1
#'
#' # Compute updated log marginal likelihood
#' log_marginal_likelihood <- update_lmm_residual_variance_cat(
#'   residual_variance, random_effect_variance, d_seq, d_seq_star, q, q_star, tau
#' )
#' log_marginal_likelihood
update_lmm_variance <- function(residual_variance,
                                random_effect_variance,
                                obs_z_eigenvalues,
                                syn_z_eigenvalues,
                                obs_adjusted_residuals,
                                syn_adjusted_residuals,
                                tau) {
  obs_log_likelihood <- sum(
    log(random_effect_variance * obs_z_eigenvalues + residual_variance) +
      obs_adjusted_residuals^2 / (random_effect_variance * obs_z_eigenvalues + residual_variance)
  )

  syn_log_likelihood <- tau / length(syn_z_eigenvalues) * sum(
    log(random_effect_variance * syn_z_eigenvalues + residual_variance) +
      syn_adjusted_residuals^2 / (random_effect_variance * syn_z_eigenvalues + residual_variance)
  )

  return(obs_log_likelihood + syn_log_likelihood)
}

validate_lmm_input <- function(
    formula,
    cat_init,
    tau = NULL,
    residual_variance_0 = NULL,
    random_effect_variance_0 = NULL,
    coefs_0 = NULL,
    optimize_domain = NULL,
    max_iter = NULL,
    tol = NULL,
    tau_seq = NULL,
    cross_validation_fold_num = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau, incl_0 = TRUE)
  validate_positive("tau_seq", tau_seq, incl_0 = TRUE, is_vector = TRUE)
  validate_positive("residual_variance_0", residual_variance_0)
  validate_positive("random_effect_variance_0", random_effect_variance_0)
  validate_positive("optimize_domain", optimize_domain, incl_0 = TRUE, is_vector = TRUE)
  validate_positive("max_iter", max_iter)
  validate_positive("tol", tol)
  validate_positive("cross_validation_fold_num", cross_validation_fold_num)

  # Check if cat_init is generated from cat_lmm_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_lmm_initialization")) {
    stop(
      patse0(
        "The provided `cat_init` must be generated from `cat_lmm_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if length of coefs_0 equal to number of columns of observation x
  if (length(coefs_0) != ncol(cat_init$obs_x)) {
    stop(
      paste0(
        "The provided `coefs_0` must have same length as the number of columns of `cat_init$obs_x`. <",
        paste(coefs_0, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if optimize_domain is length of 2
  if (length(optimize_domain) != 2) {
    stop(
      paste0(
        "The provided `optimize_domain` must have only length 2. <",
        paste(optimize_domain, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }
}

validate_lmm_initialization_input <- function(
    formula,
    data,
    X_cols,
    Y_col,
    Z_cols,
    group_col,
    syn_size) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)

  # Check if formula is for GLMs or LMM
  if (grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for Linear Mixed Model (LMM) or Generalized Linear Model (GLMs)."
      ),
      call. = FALSE
    )
  }

  # Check if the formula is overly complex
  if (length(attr(terms(formula), "term.labels")) >= (ncol(data) - 1)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if all columns has been identified and inlcuded in `data`
  colnames_subset <- c(X_cols, Y_col, Z_cols, group_col)
  if ((!all(colnames(data) %in% colnames_subset)
  ) | (length(unique(colnames_subset)) != length(colnames_subset))) {
    stop(
      paste0(
        "Please enter `X_cols`, `Y_col`, `Z_cols` and `group_col` with no overlap,",
        "and inlcude all column names showed in `data`.\n",
        "Or, please enter right group column name in the formula"
      ),
      call. = FALSE
    )
  }

  if (!is.null(syn_size) & syn_size < (4 * ncol(data))) {
    warning(
      "The provided `syn_size` should is too small, recommand to choose a larger value.",
      call. = FALSE
    )
  }
}
# ----- General -----

validate_positive <- function(param_name,
                              param_value,
                              incl_0 = FALSE,
                              is_vector = FALSE) {
  if (!is.null(param_value)) {
    condition <- ifelse(incl_0, any(param_value <= 0), any(param_value < 0))
    stop_add_on <- ifelse(incl_0, "non-negative", "positive")
    if (is_vector) {
      if (condition | length(param_value) < 1) {
        stop(
          paste0(
            "The provided `",
            param_name,
            "` should be a vector of ",
            stop_add_on,
            " numbers. <",
            paste(param_value, collapse = ", "),
            "> is invalid."
          ),
          call. = FALSE
        )
      }
    } else if (condition | length(param_value) != 1) {
      stop(
        paste0(
          "The provided `",
          param_name,
          "` should be a ",
          stop_add_on,
          " number. <",
          param_value,
          "> is invalid."
        ),
        call. = FALSE
      )
    }
  }
}

#' Get a Stan model
#'
#' This function loads and returns a Stan model based on specified parameters.
#'
#' @param type A character string specifying the model type. Choices are "glm" or "cox".
#' @param glm_family Optional character string specifying the GLM family (e.g., "binomial").
#' @param joint_tau Logical indicating if joint tau is used.
#' @param glm_binomial_joint_theta Logical indicating if joint theta is used for binomial GLM.
#' @param glm_binomial_joint_alpha Logical indicating if joint alpha is used for binomial GLM.
#' @param glm_binomial_gibbs Logical indicating if Gibbs sampling is used for binomial GLM.
#'
#' @return A Stan model object.
#'
#' @export
get_stan_model <- memoise::memoise(function(type = c("glm", "cox"),
                                            glm_family_string = NULL,
                                            joint_tau = FALSE,
                                            glm_binomial_joint_theta = FALSE,
                                            glm_binomial_joint_alpha = FALSE,
                                            glm_binomial_gibbs = FALSE,
                                            glm_gaussian_known_variance = FALSE) {
  # Ensure type is valid
  type <- match.arg(type)

  # Construct the filename based on parameters
  filename <- paste0(
    "~/catalytic/R/rstan_program/",
    paste0(
      c(
        type,
        glm_family_string,
        if (any(joint_tau, glm_binomial_joint_theta, glm_binomial_joint_alpha)) "joint",
        if (type == "glm" && glm_family_string == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs) "theta",
        if (type == "glm" && glm_family_string == "binomial" && glm_binomial_joint_alpha && !glm_binomial_gibbs) "alpha",
        if (type == "glm" && glm_family_string == "gaussian" && glm_gaussian_known_variance) "known_variance"
      ),
      collapse = "_"
    ),
    ".stan"
  )

  # Load and return the Stan model
  return(rstan::stan_model(file = filename))
})

#' Compute Linear Predictor
#'
#' This function computes the linear predictor from a matrix of predictor variables and a vector of coefficients.
#' It handles cases with and without an intercept term.
#'
#' @param X A matrix of predictor variables.
#' @param coefs A vector of coefficients. It should be either the same length as the number of columns in X
#'               (for models without an intercept) or one more than the number of columns in X (for models with an intercept).
#'
#' @return A vector of linear predictor values.
#'
#' @export
#'
#' @examples
#' X <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
#' coefs <- c(0.5, 1.5)
#' get_linear_predictor(X, coefs)
#'
#' coefs_with_intercept <- c(1, 0.5, 1.5)
#' get_linear_predictor(X, coefs_with_intercept)
get_linear_predictor <- function(X, coefs) {
  X <- as.matrix(X)
  coefs <- as.vector(coefs)

  # Case 1: Coefficients vector length matches number of columns in X (no intercept)
  if (length(coefs) == ncol(X)) {
    return(X %*% coefs)

    # Case 2: Coefficients vector length is one more than number of columns in X (with intercept)
  } else if (length(coefs) == (ncol(X) + 1)) {
    return(coefs[1] + X %*% coefs[-1])
  }

  # Error case: dimensions do not match
  stop(
    "\nCoefficient and predictor dimensions do not match:\n",
    "Coefficients:\n", coefs,
    "\nCovariance:\n", X
  )
}

#' Compute Discrepancy Measures
#'
#' This function computes various discrepancy measures between observed and estimated values.
#' It supports different methods including logarithmic error, square error, classification error, and logistic deviance.
#'
#' @param discrepancy_method A character string specifying the discrepancy method to use. Options are:
#'   \describe{
#'     \item{"logarithmic_error"}{Logarithmic error, suitable for probabilities.}
#'     \item{"square_error"}{Mean squared error.}
#'     \item{"classification_error"}{Classification error, suitable for binary outcomes.}
#'     \item{"logistic_deviance"}{Logistic deviance, computed using a GLM model.}
#'   }
#' @param family A GLM family object (e.g., "binomial") used to compute logistic deviance.
#' @param X A matrix of predictor variables.
#' @param Y A vector or data frame of observed values.
#' @param coefs A vector of coefficients for the GLM model.
#' @param est_Y A vector of estimated values. If not provided, it will be computed using `get_glm_mean` with the specified `family`.
#'
#' @return A numeric value representing the discrepancy between observed and estimated values.
#'
#' @export
#'
#' @examples
#' # Example with square error
#' X <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
#' coefs <- c(0.5, 1.5)
#' Y <- c(1, 0)
#' est_Y <- get_glm_mean("binomial", X, coefs)
#' discrepancy <- get_discrepancy("square_error", X = X, Y = Y, coefs = coefs, est_Y = est_Y)
#' discrepancy
get_discrepancy <- function(
    discrepancy_method = c(
      "logarithmic_error",
      "square_error",
      "classification_error",
      "logistic_deviance"
    ),
    family_string = NULL,
    X = NULL,
    Y = NULL,
    coefs = NULL,
    est_Y = NULL) {
  # Match the provided discrepancy method to one of the allowed methods
  discrepancy_method <- match.arg(discrepancy_method)

  # Compute estimated values (est_Y) if not provided, using get_glm_mean if family is specified
  if (is.null(est_Y) & !is.null(family_string)) {
    est_Y <- get_glm_mean(
      family_string = family_string,
      X = X,
      coefs = coefs
    )
  }

  # Apply specific preprocessing for logarithmic error
  if (discrepancy_method == "logarithmic_error") {
    Y <- pmax(0.0001, pmin(0.9999, Y))
    est_Y <- pmax(0.0001, pmin(0.9999, est_Y))
  }

  # Compute and return the discrepancy measure based on the chosen method
  return(
    switch(discrepancy_method,
      "logarithmic_error" = mean(Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y))),
      "square_error" = mean((Y - est_Y)^2),
      "classification_error" = mean(Y - 2 * Y * est_Y + est_Y),
      "logistic_deviance" = -get_glm_log_density(family_string = "binomial", X = X, Y = Y, coefs = coefs)
    )
  )
}



#' Get Adjusted Data
#'
#' This function adjusts the data in the `cat_init` object according to the specified formula.
#' It transforms the variables in the formula to numeric and separates the data into observed
#' and synthetic parts.
#'
#' @param cat_init A list containing the initial data with components:
#'   - `x`: A matrix or data frame of raw data.
#'   - `obs_size`: Number of observed samples.
#'   - `size`: Total number of samples (observed + synthetic).
#' @param formula A formula specifying the variables to be adjusted.
#'
#' @return The updated `cat_init` list with adjusted data:
#'   - `adj_x`: Adjusted data frame for all samples.
#'   - `adj_obs_x`: Adjusted data frame for observed samples.
#'   - `adj_syn_x`: Adjusted data frame for synthetic samples.
#'
#' @export
#'
#' @examples
#' cat_init <- list(
#'   x = data.frame(a = 1:5, b = 6:10),
#'   obs_size = 3,
#'   size = 5
#' )
#' formula <- "x ~ b"
#' adjusted_data <- get_adjusted_data(cat_init, formula)
#' adjusted_data
get_adjusted_cat_init <- function(cat_init,
                                  formula_rhs,
                                  type = c("glm", "cox")) {
  # Generate model frame for predictors using the right-hand side of the formula
  cat_init$adj_x <- model.frame(
    formula_rhs,
    as.data.frame(cat_init$x)
  )

  # Split adjusted data into observed and synthetic parts
  cat_init$adj_obs_x <- cat_init$adj_x[1:cat_init$obs_size, , drop = FALSE]
  cat_init$adj_syn_x <- cat_init$adj_x[(cat_init$obs_size + 1):cat_init$size, , drop = FALSE]

  return(cat_init)
}

#' Resample Data Columns and Apply Processing Methods
#'
#' This function resamples data by columns, applies specific processing methods
#' depending on the data type and characteristics, and records the resampling process.
#'
#' @param data A data frame to be resampled.
#' @param resample_size An integer specifying the size of the resample.
#' @param data_degree A numeric vector indicating the degree of each column in the data (optional).
#' @param resample_only A logical value indicating whether to return only the resampled data (default is FALSE).
#'
#' @return A list containing:
#'   \item{resampled_data}{A data frame of resampled data.}
#'   \item{resampled_data_record}{A data frame recording the resampling process for each column.}
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 3, 4, 5),
#'   x2 = factor(c("A", "B", "A", "B", "A")),
#'   x3 = c(1.1, 2.2, 3.3, 4.4, 5.5)
#' )
#' resampled <- get_resampled_data(data, resample_size = 10)
#' print(resampled$resampled_data)
#' print(resampled$resampled_data_record)
get_resampled_data <- function(data,
                               resample_size,
                               data_degree = NULL,
                               resample_only = FALSE) {
  # Resample a column
  resample_column <- function(col, size) {
    sample(col, size = size, replace = TRUE)
  }

  # Calculate skewness
  get_skewness <- function(lst) {
    n <- length(lst)
    (sum((lst - mean(lst))^3) / n) / (sum((lst - mean(lst))^2) / n)^(3 / 2)
  }

  # Resampling methods
  get_resampling_methods <- function(x,
                                     original_x,
                                     method = c(
                                       "symmetrizing",
                                       "smoothing",
                                       "flattening",
                                       "deskewing"
                                     )) {
    # Interquartile range (IQR) calculation
    get_iqr <- function(x, invl = 0.5) {
      right_quantile <- (2 - invl) / 2
      iqr_value <- 0
      while (iqr_value == 0 & right_quantile < 1) {
        iqr_value <- diff(
          quantile(
            as.numeric(x),
            c(1 - right_quantile, right_quantile),
            na.rm = FALSE,
            names = FALSE,
            type = 7
          )
        )
        right_quantile <- right_quantile + 0.05
      }

      return(iqr_value)
    }

    # Apply the selected resampling method
    switch(match.arg(method),
      "symmetrizing" = mean(original_x) + sample(
        c(-1, 1),
        length(x),
        replace = TRUE
      ) * (x - mean(original_x)),
      "smoothing" = x + stats::rnorm(
        length(x),
        0,
        min(diff(sort(unique(original_x)))) / 10
      ),
      "flattening" = replace(
        x,
        sample(length(x), length(x) / 2),
        sample(unique(original_x),
          length(x) / 2,
          replace = TRUE
        )
      ),
      "deskewing" = replace(
        x,
        sample(length(x), length(x) / 2),
        truncnorm::rtruncnorm(length(x) / 2,
          a = min(original_x),
          b = max(original_x),
          mean(original_x, trim = 0.1),
          -get_iqr(original_x) / (2 * qnorm(1 / 4))
        )
      )
    )
  }

  # Process each column based on resampling methods
  process_column <- function(column, original_column, degree) {
    process <- "Coordinate"
    data_type <- ifelse(is.continuous(original_column), "Continuous", "Category")

    if (data_type == "Continuous") {
      if (degree >= (length(unique(original_column)) - 1)) {
        process <- c(process, "Smoothing")
        column <- get_resampling_methods(
          column,
          original_column,
          "smoothing"
        )
      }

      if (abs(get_skewness(original_column)) >= 1) {
        process <- c(process, "Deskewing")
        column <- get_resampling_methods(
          column,
          original_column,
          "deskewing"
        )
      }
    } else {
      if (any(table(original_column) < length(original_column) / (2 * length(table(original_column))))) {
        process <- c(process, "Flattening")
        column <- get_resampling_methods(
          column,
          original_column,
          "flattening"
        )
      }
    }

    return(list(
      column = column,
      data_type = data_type,
      process = paste(process, collapse = "->")
    ))
  }

  # Resample the data by columns
  resampled_data <- as.data.frame(lapply(data, resample_column, size = resample_size))

  # Return resampled data only if resample_only is TRUE
  if (resample_only) {
    return(list(
      resampled_data = resampled_data,
      resampled_data_record = NULL
    ))
  }

  # Record resampling processes for each column
  resampled_data_record <- do.call(
    rbind,
    lapply(
      1:ncol(resampled_data),
      function(n) {
        processed <- process_column(resampled_data[, n], data[, n], data_degree[n])
        resampled_data[, n] <<- processed$column
        data.frame(
          Covariate = colnames(resampled_data)[n],
          Type = processed$data_type,
          Process = processed$process
        )
      }
    )
  )

  return(list(
    resampled_data = resampled_data,
    resampled_data_record = resampled_data_record
  ))
}

#' Convert Formula to String
#'
#' This function converts a formula object to a character string. It removes extra whitespace and formats the formula as a single line.
#'
#' @param formula A formula object to be converted to a string.
#'
#' @return A character string representing the formula.
#'
#' @export
#'
#' @examples
#' formula <- y ~ x1 + x2
#' get_formula_string(formula)
#'
#' formula_with_interaction <- y ~ x1 * x2
#' get_formula_string(formula_with_interaction)
get_formula_string <- function(formula) {
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  return(char)
}

#' Extract Left-Hand Side of Formula as String
#'
#' This function extracts the left-hand side (LHS) of a formula object and converts it to a character string. It uses `get_formula_string` to ensure consistent formatting.
#'
#' @param formula A formula object from which the LHS will be extracted.
#'
#' @return A character string representing the left-hand side of the formula.
#'
#' @export
#'
#' @examples
#' formula <- y ~ x1 + x2
#' get_formula_lhs(formula)
#'
#' formula_with_interaction <- y ~ x1 * x2
#' get_formula_lhs(formula_with_interaction)
get_formula_lhs <- function(formula) {
  return(get_formula_string(rlang::f_lhs(formula)))
}

#' Extract the Right-Hand Side of a Formula
#'
#' This function extracts the right-hand side (RHS) of a formula and returns it as a character string.
#' Optionally, it can include a tilde (`~`) at the beginning of the RHS.
#'
#' @param formula A formula object from which to extract the RHS.
#' @param with_tilde Logical, indicating whether to include a tilde (`~`) at the beginning of the RHS.
#'   Defaults to `FALSE`.
#'
#' @return A character string representing the right-hand side of the formula. If `with_tilde` is `TRUE`,
#'   the string includes a leading tilde.
#'
#' @export
#'
#' @examples
# formula <- y ~ x1 + x2
#' #
# get_formula_rhs(formula)
# get_formula_rhs(formula, with_tilde = TRUE)
get_formula_rhs <- function(formula, with_tilde = FALSE) {
  if (with_tilde) {
    return(paste("~", deparse(rlang::f_rhs(formula)), collapse = " "))
  }
  return(get_formula_string(rlang::f_rhs(formula)))
}

#' Check if a Variable is Continuous
#'
#' This function checks whether a given vector represents a continuous variable. A continuous variable is numeric and has more than two unique values.
#'
#' @param lst A vector to be checked.
#'
#' @return A logical value indicating whether the input vector is considered continuous. Returns `TRUE` if the vector is numeric and has more than two unique values; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with a continuous variable
#' continuous_vector <- c(1, 2, 3, 4, 5)
#' is.continuous(continuous_vector) # TRUE
#'
#' # Example with a categorical variable
#' categorical_vector <- c("A", "B", "A", "B")
#' is.continuous(categorical_vector) # FALSE
#'
#' # Example with a binary variable
#' binary_vector <- c(0, 1, 1, 0)
#' is.continuous(binary_vector) # FALSE
is.continuous <- function(lst) {
  # Check if the input 'lst' is numeric and has more than 2 unique values
  return(is.numeric(lst) & length(unique(lst)) > 2)
}

#' Standardize Data
#'
#' This function standardizes a dataset by converting columns to numeric or factor types and replacing NA values.
#' For continuous variables, NA values are replaced with either a specific numeric value or a computed statistic.
#' For categorical variables, NA values are replaced with the mode of the column.
#'
#' @param data A data frame to be standardized.
#' @param na_replace A function or numeric value used to replace NA values. If a function, it should take a vector and return a replacement value.
#'   If a numeric value, it is used directly to replace NA values in continuous columns.
#'   The default is `stats::na.omit`, which omits rows with NA values (used as an indicator here, not the actual replacement value).
#'
#' @return A data frame where columns have been converted to numeric or factor types, and NA values have been replaced according to the method specified.
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#'   numeric_col = c(1, 2, NA, 4),
#'   factor_col = factor(c("A", "B", "A", NA)),
#'   stringsAsFactors = FALSE
#' )
#' standardized_data <- get_standardized_data(data, na_replace = 0)
#' standardized_data
get_standardized_data <- function(data, na_replace = stats::na.omit) {
  # Function to convert columns to numeric or factor
  convert_column <- function(x) {
    if (is.continuous(stats::na.omit(x))) as.numeric(x) else x
  }

  # Function to replace NA values based on data type
  replace_na <- function(x) {
    if (anyNA(x)) {
      replace_value <- if (is.continuous(stats::na.omit(x))) {
        if (is.numeric(na_replace)) na_replace else na_replace(x, na.rm = TRUE)
      } else {
        unique_vals <- unique(x)
        unique_vals[which.max(tabulate(match(x, unique_vals)))]
      }

      x[is.na(x)] <- replace_value
    }
    x
  }

  # Convert columns and replace NA values
  data <- lapply(data, convert_column)
  data <- as.data.frame(if (identical(na_replace, stats::na.omit)) {
    stats::na.omit(data)
  } else {
    lapply(data, replace_na)
  })

  return(data)
}

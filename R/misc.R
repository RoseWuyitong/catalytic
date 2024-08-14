# ---- GLM Functions ----

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
    "square_error" = 2 * get_glm_mean("gaussian", X, coefs),
    "classification_error" = 2 * get_glm_mean("binomial", X, coefs),
    "logistic_deviance" = get_linear_predictor(X, coefs)
  ))
}

#' Compute Diagonal Approximate Covariance Matrix
#'
#' This function computes the diagonal elements of the approximate covariance matrix
#' for the coefficients in a generalized linear model (GLM). The covariance is derived
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
get_glm_diag_approx_cov <- function(X, model) {
  # Function to compute the second derivative matrix
  get_second_derivative <- function(X,
                                    fitted_values,
                                    prior_weights) {
    return((t(cbind(1, X))) %*%
      diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
      cbind(1, X))
  }

  # Compute the diagonal of the inverse of the Hessian matrix
  return(diag(chol2inv(chol(get_second_derivative(
    X = as.matrix(X),
    fitted_values = model$fitted.values,
    prior_weights = model$prior.weights
  )))))
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
#' family_name <- get_glm_family("gaussian")
#' family_name
#'
#' # Retrieve the name of the Poisson family with its link function
#' family_name_with_link <- get_glm_family("poisson", with_link = TRUE)
#' family_name_with_link
get_glm_family <- function(family, with_link = FALSE) {
  # Validate and retrieve the family object
  validate_family_glm <- function(f) {
    if (is.character(f)) {
      tryCatch(
        {
          f <- get(f, mode = "function", envir = parent.frame(2))
        },
        error = function(e) stop("'family' must be a valid GLM family.", call. = FALSE)
      )
    }
    if (is.function(f)) {
      f <- f()
    }
    if (!is(f, "family")) {
      stop("'family' must be a valid GLM family.", call. = FALSE)
    }
    return(f)
  }

  # Retrieve and validate the family object
  family <- validate_family_glm(family)

  # Return the family name or family name with link function
  if (with_link) {
    return(paste0(family$family, " [", family$link, "]"))
  } else {
    return(family$family)
  }
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
get_glm_mean <- function(family, X, coefs) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X, coefs)

  # Return the mean based on the specified family
  return(switch(get_glm_family(family),
    "binomial" = 1 / (1 + exp(-lp)),
    "gaussian" = lp
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
get_glm_log_density <- function(family, X, Y, coefs, weights = 1) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X, coefs)

  # Ensure Y is a numeric vector
  Y <- if (is.data.frame(Y)) {
    as.vector(apply(Y, 2, function(x) as.numeric(as.character(x))))
  } else {
    as.vector(Y)
  }

  # Calculate the log density based on the specified family
  return(
    switch(get_glm_family(family),
      "binomial" = sum(weights * (Y * lp - log(1 + exp(lp)))),
      "gaussian" = 0.5 * sum(weights * (log(2 * pi * var(Y - lp)) + 1))
    )
  )
}

#' Compute Gradient of Log Density for GLM Families
#'
#' This function calculates the gradient of the log density with respect to the coefficients
#' for a given GLM family based on the provided predictors, response variable, and weights.
#'
#' @param family Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param Y Vector. The response variable.
#' @param coefs Numeric vector. The coefficients for the GLM.
#' @param weights Numeric vector. The weights for the GLM. Default is 1.
#'
#' @return Numeric vector. The gradient of the log density with respect to the coefficients.
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
get_glm_log_density_grad <- function(family,
                                     X,
                                     Y,
                                     coefs,
                                     weights = 1) {
  # Compute the mean based on the GLM family
  mean_val <- get_glm_mean(family, X, coefs)

  # Calculate the gradient of the log density
  gradient <- weights * (as.vector(Y) - mean_val)
  return(c(t(gradient) %*% as.matrix(cbind(1, X))))
}

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
get_glm_sample_data <- function(family = "gaussian",
                                n = 10,
                                mean = 0,
                                sd = 1) {
  switch(get_glm_family(family),
    "binomial" = stats::rbinom(n, 1, mean),
    "gaussian" = stats::rnorm(n, mean, sd)
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
get_glm_bayes_joint_suggestion <- function(family,
                                           stan_iter,
                                           alpha,
                                           stan_sample_model,
                                           binomial_use_theta,
                                           binomial_joint_alpha,
                                           binomial_tau_lower) {
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

  # Check the GLM family and proceed with specific handling for binomial
  if (get_glm_family(family) == "binomial") {
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
        } else if (alpha < 5) {
          cat("\nSUGGESTION: Please consider using a larger `alpha` (e.g.,", min(alpha + 0.1, 5), ")\n")
        } else {
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
}

#' Check Dataset Size and Risk Estimate Method
#'
#' This function checks the size of the observational dataset against
#' the specified risk estimate method. It provides warnings or prompts
#' to guide users in selecting appropriate risk estimation methods based
#' on the dataset size.
#'
#' @param obs_size Numeric. The size of the observational dataset.
#' @param risk_estimate_method Character. The method used for risk estimation.
#'        Options are `"parametric_bootstrap"`, `"steinian_estimate"`, or `"cross_validation"`.
#'
#' @return NULL. The function stops or issues warnings based on the dataset size
#'         and the chosen risk estimate method.
#'
#' @export
check_glm_tune_size_and_risk_estimate_method <- function(obs_size,
                                                         risk_estimate_method) {
  # Check if the dataset size exceeds the threshold for the given risk estimate method
  if (obs_size > if (risk_estimate_method == "parametric_bootstrap") {
    600
  } else if (risk_estimate_method == "steinian_estimate") {
    400
  } else {
    Inf
  }) {
    # Prompt the user to continue if the dataset size is large
    response <- menu(c("YES", "NO"), title = paste0(
      "The large size of the dataset ",
      "could make estimate time longer than expected.",
      "\nDo you want to continue?\n"
    ))

    # Stop the process if the user chooses "NO"
    if (response == 2) {
      stop("Process terminated by user.", call. = FALSE)
    }
  }

  # Issue a warning if the dataset size is larger than recommended for the risk estimate method
  if (obs_size > if (risk_estimate_method == "parametric_bootstrap") {
    300
  } else if (risk_estimate_method == "steinian_estimate") {
    200
  } else {
    Inf
  }) {
    warning(paste(
      "Due to the large size of the dataset,",
      "Please consider utilizing `risk_estimate_method = cross_validation`",
      "for quicker estimate process."
    ), call. = FALSE)
  }

  # Issue a warning if the dataset size is insufficient for cross-validation
  if (obs_size < if (risk_estimate_method == "cross_validation") 30 else 0) {
    warning(paste0(
      "The size of the observational dataset is insufficient ",
      "for cross validation to yield optimal results.",
      "\nSUGGESTION: Please consider utilizing alternative methods",
      "(e.g. risk_estimate_method = parametric_bootstrap) ",
      "for better estimate."
    ), call. = FALSE)
  }
}


# ---- COX Functions ----


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
get_cox_qr_solve <- function(coef_matrix, constant_vector) {
  # Calculate solution using QR decomposition with error handling
  qr_solution <- tryCatch(
    {
      # Attempt to solve the linear system
      solve(coef_matrix, constant_vector)
    },
    error = function(e) {
      # In case of error, return a vector of default values
      rep(0.0001, length(constant_vector))
    }
  )

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
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
                            entry_points = NULL) {
  # Convert inputs to appropriate formats
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Default entry points to a vector of zeros if not provided
  if (is.null(entry_points)) {
    entry_points <- rep(0, nrow(X))
  }

  # Identify the indices of events (status == 1)
  censored_status <- which(status == 1)

  # Initialize the Hessian matrix
  p <- ncol(X)
  hessian <- matrix(0, p, p)

  # Calculate the exponential of the linear predictor
  exp_risk_lp <- exp(get_linear_predictor(X, coefs))

  # Compute the Hessian matrix contributions for each event
  for (i in censored_status) {
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
    risk_set <- X[risk_set_idx, , drop = FALSE]
    exp_risk_lp_risk_set <- exp_risk_lp[risk_set_idx]
    sum_exp_risk_lp <- sum(exp_risk_lp_risk_set)

    # Compute the Hessian matrix contribution
    risk_exp_risk_lp <- t(risk_set) %*% exp_risk_lp_risk_set
    risk_sqrt_exp_risk_lp <- t(sweep(risk_set, 1, sqrt(exp_risk_lp_risk_set), "*"))

    hessian <- hessian +
      tcrossprod(risk_sqrt_exp_risk_lp) / sum_exp_risk_lp -
      tcrossprod(risk_exp_risk_lp) / sum_exp_risk_lp^2
  }

  # Return the negative Hessian
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
get_cox_syn_hessian <- function(X, time, coefs, hazard_constant) {
  # Convert inputs to appropriate formats
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate the linear predictors
  syn_lp <- as.numeric(get_linear_predictor(X, coefs))

  # Compute the Hessian using a vectorized approach
  hessian <- -hazard_constant *
    t(X) %*% (sweep(X, 1, time * exp(syn_lp), "*"))

  return(hessian)
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
                             entry_points = NULL) {
  # Convert inputs to appropriate formats
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Set entry points to zero if not provided
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  # Initialize gradient matrix
  gradient <- matrix(NA, nrow = ncol(X), ncol = nrow(X))

  # Loop over each observation to compute the gradient
  for (i in seq_len(nrow(X))) {
    # Get the indices of the risk set for the current observation
    risk_set_idx <- get_cox_risk_set_idx(time[i], entry_points, time, status)

    # Skip if the risk set is empty
    if (length(risk_set_idx) == 0) next

    # Extract risk set covariates and compute exponential risk
    risk_set <- as.matrix(X[risk_set_idx, ])
    exp_risk_set <- exp(risk_set %*% coefs)

    # Compute the gradient for the current observation
    gradient[, i] <- X[i, ] -
      (t(risk_set) %*% exp_risk_set / sum(exp_risk_set))
  }

  # Return the gradient vector by summing over all observations weighted by the status vector
  return(as.vector(gradient %*% status))
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
  # Convert inputs to appropriate formats
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate linear predictors
  syn_lp <- get_linear_predictor(X, coefs)

  # Calculate the expression exp(syn_lp) once for all rows
  exp_syn_lp <- exp(syn_lp)

  # Calculate the gradient in a vectorized manner
  gradient <- colSums((1 - hazard_constant * time * exp_syn_lp) * X)

  return(as.numeric(gradient))
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
                                       entry_points = NULL) {
  # Convert inputs to appropriate formats
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  # Find indices of censored observations
  censored_status <- which(as.numeric(status) == 1)

  # Initialize partial likelihood
  pl <- 0

  # Calculate partial likelihood for each censored observation
  for (i in censored_status) {
    # Get indices of the risk set for the current observation
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    # Compute the linear predictors for the risk set
    exp_risk_lp <- exp(as.vector(get_linear_predictor(
      matrix(as.matrix(X[risk_set_idx, ]), nrow = length(risk_set_idx)),
      coefs
    )))

    # Update partial likelihood
    pl <- pl + sum(X[i, ] * coefs) - log(sum(exp_risk_lp))
  }

  return(pl)
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
                          hazard_constant = NULL,
                          tol = 1e-5,
                          max_iter = 25,
                          entry_points = NULL,
                          initial_values = NULL) {
  # Ensure data is in the correct format
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Calculate hazard_constant if not provided
  if (is.null(hazard_constant)) hazard_constant <- sum(status) / mean(time) / nrow(X)

  # Set default entry points if not provided
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  # Initialize coefficients
  coefs <- if (is.null(initial_values)) rep(0, ncol(X)) else initial_values

  iter <- 1
  diff_new_old <- TRUE

  # Iteratively estimate coefficients
  while ((iter < max_iter) & diff_new_old) {
    coefs_old <- coefs

    # Calculate linear predictor
    exp_lp <- exp(get_linear_predictor(X, coefs_old))

    # Calculate gradient and Hessian
    gd_old <- get_cox_syn_gradient(X, time, coefs_old, hazard_constant)
    hess_old <- get_cox_syn_hessian(X, time, coefs_old, hazard_constant)

    # Update coefficients using QR solve
    coefs <- coefs_old - get_cox_qr_solve(hess_old, gd_old)

    iter <- iter + 1

    # Stop if NA values are encountered
    if (any(is.na(coefs))) stop("Please note NA appears in iteration.")

    # Stop if coefficients become too large
    if (norm(as.numeric(coefs), "2") > 1000000) stop("Please note the Euclidean norm of coefficients is too large.")

    # Check for convergence
    diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > tol * norm(as.numeric(coefs_old), "2")
  }

  # Calculate the final linear predictor and return the kappa value
  lp_new <- as.numeric(get_linear_predictor(X, coefs))
  return(sum(lp_new - exp(lp_new) * hazard_constant * time) / nrow(X))
}

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
  # Ensure inputs are vectors
  time_vector <- as.vector(time_vector)
  status_vector <- as.vector(status_vector)
  indicator_vector <- as.vector(indicator_vector)

  # Function to find the index where an event occurs or reaches censoring time
  find_occur_idx <- function(time) {
    min(which(indicator_vector >= time))
  }

  # Identify the indices of events (status == 1)
  censored_status <- which(status_vector == 1)

  # Initialize risk and failure sets as zero matrices
  risk_set <- failure_set <- matrix(
    0,
    length(time_vector),
    length(indicator_vector)
  )

  # Process subjects who experienced an event (status = 1)
  for (i in censored_status) {
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
  # Convert inputs to appropriate vector formats
  time_vector <- as.vector(time_vector)
  entry_vector <- as.vector(entry_vector)
  status_vector <- as.vector(status_vector)

  # Find indices where subjects are at risk at the given time of interest
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & status_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}


# ---- LMM Functions ----

#' Update the log marginal likelihood for sigma^2 in a hierarchical linear mixed model with catalytic prior
#'
#' This function calculates the log marginal likelihood for sigma^2 in a hierarchical linear mixed model (LMM)
#' with a catalytic prior. The marginal likelihood is updated based on provided parameters and data sequences.
#'
#' @param sigma_square A scalar representing the variance of the residual errors.
#' @param sigma_u_square A scalar representing the variance of the random effects.
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
#' sigma_square <- 1
#' sigma_u_square <- 0.5
#' d_seq <- c(2, 3)
#' d_seq_star <- c(1, 2)
#' q <- c(0.5, 1)
#' q_star <- c(0.3, 0.7)
#' tau <- 0.1
#'
#' # Compute updated log marginal likelihood
#' log_marginal_likelihood <- update_lmm_sigma_square_cat(
#'   sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau
#' )
#' log_marginal_likelihood
update_lmm_sigma_square_cat <- function(sigma_square,
                                        sigma_u_square,
                                        d_seq,
                                        d_seq_star,
                                        q,
                                        q_star,
                                        tau) {
  # Ensure d_seq and d_seq_star are vectors, and q and q_star are vectors of the same length as d_seq and d_seq_star respectively
  d_seq <- as.numeric(d_seq)
  d_seq_star <- as.numeric(d_seq_star)
  q <- as.numeric(q)
  q_star <- as.numeric(q_star)

  # Compute the marginal likelihood for the data and prior
  log_likelihood_data <- sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
  log_likelihood_prior <- tau / length(d_seq_star) * sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))

  return(log_likelihood_data + log_likelihood_prior)
}

# ---- Other Functions ----

#' Compute the H matrix for the catalytic prior
#'
#' This function calculates the H matrix for the catalytic prior given observed and synthetic data matrices,
#' and a precision parameter \code{tau}. The H matrix is used in the context of Bayesian analysis with a
#' catalytic prior distribution.
#'
#' @param obs_x A matrix of observed data (size: n_observations x n_features).
#' @param syn_x A matrix of synthetic data (size: n_synthetic_observations x n_features).
#' @param tau A scalar representing the precision parameter of the catalytic prior.
#'
#' @return A matrix representing the computed H matrix.
#'
#' @export
#'
#' @examples
#' obs_x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' syn_x <- matrix(c(7, 8, 9, 10, 11, 12), ncol = 2)
#' tau <- 0.1
#'
#' # Compute the H matrix
#' H_matrix <- get_h_matrix(obs_x, syn_x, tau)
#' H_matrix
get_h_matrix <- function(obs_x, syn_x, tau) {
  obs_x <- as.matrix(obs_x)
  syn_x <- as.matrix(syn_x)
  obs_size <- nrow(obs_x)
  syn_size <- nrow(syn_x)

  # Compute the H matrix
  return(obs_x %*% (
    chol2inv(chol((t(obs_x) %*% obs_x + (tau / syn_size) * (t(syn_x) %*% syn_x))))) %*% (
    t(obs_x) + tau / (syn_size * obs_size) * t(syn_x) %*% rep(1, syn_size) %*% t(rep(1, obs_size))
  ))
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
                                            glm_family = NULL,
                                            joint_tau = FALSE,
                                            glm_binomial_joint_theta = FALSE,
                                            glm_binomial_joint_alpha = FALSE,
                                            glm_binomial_gibbs = FALSE) {
  # Ensure type is valid
  type <- match.arg(type)

  # Get the GLM family if provided
  glm_family <- if (!is.null(glm_family)) get_glm_family(glm_family) else NULL

  # Construct the filename based on parameters
  filename <- paste0(
    "~/catalytic/R/rstan_script/",
    paste0(
      c(
        type,
        glm_family,
        if (any(joint_tau, glm_binomial_joint_theta, glm_binomial_joint_alpha)) "joint" else NULL,
        if (type == "glm" && glm_family == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs) "theta" else NULL,
        if (type == "glm" && glm_family == "binomial" && glm_binomial_joint_alpha && !glm_binomial_gibbs) "alpha" else NULL
      ),
      collapse = "_"
    ),
    ".stan"
  )

  # Load and return the Stan model
  return(rstan::stan_model(file = filename))
})


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
get_adjusted_data <- function(cat_init, formula) {
  # Convert the formula to a model frame using the right-hand side
  cat_init$adj_x <- model.frame(
    get_formula_rhs(as.formula(formula), TRUE),
    as.data.frame(cat_init$x)
  )

  # Ensure all variables are numeric
  cat_init$adj_x <- as.data.frame(apply(cat_init$adj_x, 2, function(x) as.numeric(as.character(x))))

  # Preserve the row names from the original data
  rownames(cat_init$adj_x) <- rownames(cat_init$x)

  # Split the adjusted data into observed and synthetic parts
  cat_init$adj_obs_x <- as.data.frame(cat_init$adj_x[1:cat_init$obs_size, ])
  cat_init$adj_syn_x <- as.data.frame(cat_init$adj_x[(cat_init$obs_size + 1):cat_init$size, ])

  # Ensure column names are consistent
  colnames(cat_init$adj_obs_x) <- colnames(cat_init$adj_syn_x) <- colnames(cat_init$adj_x)

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
  # Function to resample a column
  resample_column <- function(column, size) {
    if (is.factor(column)) {
      factor(sample(column, size = size, replace = TRUE), levels = levels(column))
    } else if (is.numeric(column) || is.integer(column)) {
      as.numeric(sample(column, size = size, replace = TRUE))
    } else {
      sample(column, size = size, replace = TRUE)
    }
  }

  # Function to calculate skewness
  get_skewness <- function(lst) {
    n <- length(lst)
    return((sum((lst - mean(lst))^3) / n) / (sum((lst - mean(lst))^2) / n)^(3 / 2))
  }

  # Function to determine the resampling process for each column
  process_column <- function(column, original_column, degree) {
    process <- "Coordinate"
    data_type <- if (is.continuous(original_column)) "Continuous" else "Category"

    if (data_type == "Continuous") {
      if (degree >= (length(unique(original_column)) - 1)) {
        process <- c(process, "Smoothing")
        column <- do.call(get_resampling_methods, list(column, original_column, "smoothing"))
      }
      if (abs(get_skewness(original_column)) >= 1) {
        process <- c(process, "Deskewing")
        column <- do.call(get_resampling_methods, list(column, original_column, "deskewing"))
      }
    } else {
      tabl <- table(original_column)
      if (any(tabl < length(original_column) / (2 * length(tabl)))) {
        process <- c(process, "Flattening")
        column <- do.call(get_resampling_methods, list(column, original_column, "flattening"))
      }
    }

    return(list(
      column = column,
      data_type = data_type,
      process = paste(process, collapse = "->")
    ))
  }


  data <- as.data.frame(data)
  if (is.null(data_degree)) data_degree <- rep(1, ncol(data))

  # Resample the data by columns
  resampled_data <- as.data.frame(lapply(data, resample_column, size = resample_size))

  # Return resampled data only if resample_only is TRUE
  if (resample_only) {
    return(list(resampled_data = resampled_data, resampled_data_record = NULL))
  }

  # Process each column to determine resampling methods applied
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

#' Get Resampling Methods
#'
#' This function applies various resampling methods to a given dataset to transform it.
#' The methods include symmetrizing, smoothing, flattening, and deskewing.
#'
#' @param x A numeric vector to be resampled.
#' @param original_x The original numeric vector from which `x` was derived.
#' @param method A character vector specifying the resampling method. Options are:
#'   - `"symmetrizing"`: Symmetrizes `x` around the mean of `original_x`.
#'   - `"smoothing"`: Adds normal noise to `x` based on the unique values of `original_x`.
#'   - `"flattening"`: Replaces half of `x` with random values from `original_x`.
#'   - `"deskewing"`: Adjusts `x` to reduce skewness using a truncated normal distribution.
#'
#' @return A numeric vector of the same length as `x` after applying the specified resampling method.
#'
#' @export
#'
#' @examples
#' original_x <- rnorm(10)
#' x <- original_x + rnorm(10, 0, 0.5)
#' resampled_x <- get_resampling_methods(x, original_x, method = "smoothing")
#' resampled_x
get_resampling_methods <- function(x,
                                   original_x,
                                   method = c(
                                     "symmetrizing",
                                     "smoothing",
                                     "flattening",
                                     "deskewing"
                                   )) {
  # Function to calculate the interquartile range (IQR) with an adjustable interval
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

  # Match the resampling method
  method <- match.arg(method)

  # Apply the selected resampling method
  switch(method,
    "symmetrizing" = mean(original_x) +
      sample(c(-1, 1), length(x), replace = TRUE) *
        (x - mean(original_x)),
    "smoothing" = x +
      stats::rnorm(length(x), 0, min(diff(sort(unique(original_x)))) / 10),
    "flattening" = {
      x[sample(length(x), length(x) / 2)] <- sample(
        unique(original_x),
        length(x) / 2,
        replace = TRUE
      )
      return(x)
    },
    "deskewing" = {
      x[sample(length(x), length(x) / 2)] <- truncnorm::rtruncnorm(
        length(x) / 2,
        a = min(original_x),
        b = max(original_x),
        mean(original_x, trim = 0.1),
        -get_iqr(original_x) / (2 * qnorm(1 / 4))
      )
      return(x)
    }
  )
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
  convert_to_numeric_or_factor <- function(x) {
    if (is.continuous(stats::na.omit(x))) {
      return(as.numeric(x))
    } else {
      return(x)
    }
  }

  # Function to replace NA values based on data type
  replace_df_na <- function(x) {
    get_mode <- function(lst) {
      uniq_vals <- unique(lst)
      return(uniq_vals[which.max(tabulate(match(lst, uniq_vals)))])
    }

    if (anyNA(x)) {
      if (is.continuous(stats::na.omit(x))) {
        replace_value <- if (is.numeric(na_replace)) na_replace else na_replace(x, na.rm = TRUE)
        x[is.na(x)] <- replace_value
      } else {
        x[is.na(x)] <- get_mode(x)
      }
    }
    return(x)
  }

  # Convert all columns to numeric or factor
  data <- as.data.frame(lapply(data, convert_to_numeric_or_factor))

  # Replace NA values in the data
  data <- as.data.frame(
    if (identical(stats::na.omit, na_replace)) {
      stats::na.omit(data)
    } else {
      lapply(data, replace_df_na)
    }
  )

  return(data)
}

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
  # Convert X to matrix to ensure correct multiplication
  X <- as.matrix(X)

  # Convert coefs to numeric vector
  coefs <- as.numeric(coefs)

  # Case 1: Coefficients vector length matches number of columns in X (no intercept)
  if (length(coefs) == ncol(X)) {
    # Compute linear predictor: X %*% coefs
    return(as.vector(X %*% coefs))

    # Case 2: Coefficients vector length is one more than number of columns in X (with intercept)
  } else if (length(coefs) == (ncol(X) + 1)) {
    # Compute linear predictor: intercept + X %*% coefs[-1]
    return(as.vector(coefs[1] + X %*% coefs[-1]))

    # Error case: dimensions do not match
  } else {
    stop(
      "Dimensions of coefs and predictor variables do not match.\n\n Input Coefficients:\n",
      coefs,
      "\n\nInput Covariance:\n",
      X
    )
  }
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
get_discrepancy <- function(discrepancy_method = c(
                              "logarithmic_error",
                              "square_error",
                              "classification_error",
                              "logistic_deviance"
                            ),
                            family = NULL,
                            X = NULL,
                            Y = NULL,
                            coefs = NULL,
                            est_Y = NULL) {
  # Match the provided discrepancy method to one of the allowed methods
  discrepancy_method <- match.arg(discrepancy_method)

  # Convert Y to a numeric vector if it's a data frame
  Y <- if (is.data.frame(Y)) {
    as.vector(apply(Y, 2, function(x) as.numeric(as.character(x))))
  } else {
    as.vector(Y)
  }

  # Compute estimated values (est_Y) if not provided, using get_glm_mean if family is specified
  if (is.null(est_Y) & !is.null(family)) est_Y <- get_glm_mean(family, X, coefs)

  # Apply specific preprocessing for logarithmic error
  if (discrepancy_method == "logarithmic_error") {
    # Clip Y and est_Y to avoid log(0) issues by setting a small lower bound
    Y <- pmax(0.0001, pmin(0.9999, Y))
    est_Y <- pmax(0.0001, pmin(0.9999, est_Y))
  }

  # Compute and return the discrepancy measure based on the chosen method
  return(
    switch(discrepancy_method,
      "logarithmic_error" = mean(Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y))),
      "square_error" = mean((Y - est_Y)^2),
      "classification_error" = mean(Y - 2 * Y * est_Y + est_Y),
      "logistic_deviance" = -get_glm_log_density("binomial", X, Y, coefs)
    )
  )
}


# ---- Identification Functions ----

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


#' Check if an Object is a Formula
#'
#' This function checks whether a given object is a formula by looking for the presence of the tilde (`~`) character, which is a standard indicator of formulas in R.
#'
#' @param f An object to be checked.
#'
#' @return A logical value indicating whether the input object is a formula. Returns `TRUE` if the object contains the tilde (`~`) character; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with a formula
#' formula_obj <- y ~ x + z
#' is.formula(formula_obj) # TRUE
#'
#' # Example with a non-formula object
#' non_formula_obj <- c(1, 2, 3)
#' is.formula(non_formula_obj) # FALSE
#'
#' # Example with a formula string
#' formula_string <- "y ~ x + z"
#' is.formula(formula_string) # TRUE
is.formula <- function(f) {
  return(grepl("~", deparse(f)))
}

#' Check if a Formula is for Generalized Linear Models (GLM)
#'
#' This function checks whether a given formula is suitable for Generalized Linear Models (GLM) by verifying that it contains a tilde (`~`) and is not a formula for Linear Mixed Models (LMM) or Survival Analysis.
#'
#' @param f A formula or object to be checked.
#'
#' @return A logical value indicating whether the input formula is appropriate for GLM. Returns `TRUE` if the formula contains the tilde (`~`) and is not an LMM or Survival Analysis formula; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with a GLM formula
#' glm_formula <- y ~ x + z
#' is.formula_glm(glm_formula) # TRUE
#'
#' # Example with an LMM formula
#' lmm_formula <- y ~ x + (1 | group)
#' is.formula_glm(lmm_formula) # FALSE
#'
#' # Example with a Survival Analysis formula
#' surv_formula <- Surv(time, status) ~ x
#' is.formula_glm(surv_formula) # FALSE
#'
#' # Example with a non-formula object
#' non_formula_obj <- c(1, 2, 3)
#' is.formula_glm(non_formula_obj) # FALSE
is.formula_glm <- function(f) {
  if (is.formula_lmm(f) | is.formula_surv(f)) {
    return(FALSE)
  }
  return(is.formula(f))
}

#' Check if the Input is a Linear Mixed Model Formula
#'
#' This function checks if the input object is a formula suitable for linear mixed models (LMMs).
#' It returns `TRUE` if the formula contains terms with random effects, typically indicated by parentheses (`( )`);
#' otherwise, it returns `FALSE`.
#'
#' @param f An object to be checked.
#'
#' @return A logical value indicating whether the input object is a formula suitable for LMMs. Returns `TRUE` if the formula contains random effects terms (indicated by parentheses); otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with an LMM formula
#' formula_lmm <- as.formula("response ~ predictor + (1 | group)")
#' is.formula_lmm(formula_lmm) # Should return TRUE
#'
#' # Example with a GLM formula
#' formula_glm <- as.formula("response ~ predictor")
#' is.formula_lmm(formula_glm) # Should return FALSE
is.formula_lmm <- function(f) {
  # Check if the formula contains parentheses indicating random effects
  return(grepl("\\(.*\\|.*\\)", deparse(f)))
}

#' Check if the Input is a Survival Formula
#'
#' This function checks if the input object is a formula suitable for survival analysis.
#' It returns `TRUE` if the formula contains the `Surv` function, indicating it is used for survival models;
#' otherwise, it returns `FALSE`.
#'
#' @param f An object to be checked.
#'
#' @return A logical value indicating whether the input object is a formula suitable for survival analysis. Returns `TRUE` if the formula contains the `Surv` function; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with a survival formula
#' formula_surv <- as.formula("Surv(time, status) ~ predictor")
#' is.formula_surv(formula_surv) # Should return TRUE
#'
#' # Example with a GLM formula
#' formula_glm <- as.formula("response ~ predictor")
#' is.formula_surv(formula_glm) # Should return FALSE
#'
#' # Example with an LMM formula
#' formula_lmm <- as.formula("response ~ predictor + (1 | group)")
#' is.formula_surv(formula_lmm) # Should return FALSE
is.formula_surv <- function(f) {
  return(grepl("Surv", deparse(f)))
}


#' Check if Object is a Cat Initialization for GLM
#'
#' This function checks if the input object is an instance of a specific class
#' used for initializing a generalized linear model (GLM) within the catalytic framework.
#' It returns `TRUE` if the object inherits from the `"cat_initialization"` class and has
#' a `function_name` attribute equal to `"cat_glm_initialization"`.
#'
#' @param x An object to be checked.
#'
#' @return A logical value indicating whether the input object is a cat initialization for GLM.
#' Returns `TRUE` if the object inherits from the `"cat_initialization"` class and has the
#' `function_name` attribute set to `"cat_glm_initialization"`; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example of creating a cat_initialization object
#' init_obj <- list(function_name = "cat_glm_initialization")
#' class(init_obj) <- "cat_initialization"
#' is.cat_glm_init(init_obj) # Should return TRUE
#'
#' # Example with a different function_name
#' init_obj_diff <- list(function_name = "cat_other_initialization")
#' class(init_obj_diff) <- "cat_initialization"
#' is.cat_glm_init(init_obj_diff) # Should return FALSE
#'
#' # Example with a non-cat_initialization object
#' non_cat_obj <- list(function_name = "cat_glm_initialization")
#' is.cat_glm_init(non_cat_obj) # Should return FALSE
is.cat_glm_init <- function(x) {
  return(inherits(x, "cat_initialization") & x$function_name == "cat_glm_initialization")
}

#' Check if Object is a Cat Initialization for LMM
#'
#' This function checks if the input object is an instance of a specific class
#' used for initializing a linear mixed model (LMM) within the catalytic framework.
#' It returns `TRUE` if the object inherits from the `"cat_initialization"` class and has
#' a `function_name` attribute equal to `"cat_lmm_initialization"`.
#'
#' @param x An object to be checked.
#'
#' @return A logical value indicating whether the input object is a cat initialization for LMM.
#' Returns `TRUE` if the object inherits from the `"cat_initialization"` class and has the
#' `function_name` attribute set to `"cat_lmm_initialization"`; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example of creating a cat_initialization object for LMM
#' init_obj <- list(function_name = "cat_lmm_initialization")
#' class(init_obj) <- "cat_initialization"
#' is.cat_lmm_init(init_obj) # Should return TRUE
#'
#' # Example with a different function_name
#' init_obj_diff <- list(function_name = "cat_other_initialization")
#' class(init_obj_diff) <- "cat_initialization"
#' is.cat_lmm_init(init_obj_diff) # Should return FALSE
#'
#' # Example with a non-cat_initialization object
#' non_cat_obj <- list(function_name = "cat_lmm_initialization")
#' is.cat_lmm_init(non_cat_obj) # Should return FALSE
is.cat_lmm_init <- function(x) {
  return(inherits(x, "cat_initialization") & x$function_name == "cat_lmm_initialization")
}

#' Check if Object is a Cat Initialization for Cox Regression
#'
#' This function checks if the input object is an instance of a specific class
#' used for initializing a Cox regression model within the catalytic framework.
#' It returns `TRUE` if the object inherits from the `"cat_initialization"` class and has
#' a `function_name` attribute equal to `"cat_cox_initialization"`.
#'
#' @param x An object to be checked.
#'
#' @return A logical value indicating whether the input object is a cat initialization for Cox regression.
#' Returns `TRUE` if the object inherits from the `"cat_initialization"` class and has the
#' `function_name` attribute set to `"cat_cox_initialization"`; otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example of creating a cat_initialization object for Cox regression
#' init_obj <- list(function_name = "cat_cox_initialization")
#' class(init_obj) <- "cat_initialization"
#' is.cat_cox_init(init_obj) # Should return TRUE
#'
#' # Example with a different function_name
#' init_obj_diff <- list(function_name = "cat_other_initialization")
#' class(init_obj_diff) <- "cat_initialization"
#' is.cat_cox_init(init_obj_diff) # Should return FALSE
#'
#' # Example with a non-cat_initialization object
#' non_cat_obj <- list(function_name = "cat_cox_initialization")
#' is.cat_cox_init(non_cat_obj) # Should return FALSE
is.cat_cox_init <- function(x) {
  return(inherits(x, "cat_initialization") & x$function_name == "cat_cox_initialization")
}

#' Check Subset Membership
#'
#' This function checks if all elements of a `whole_set` are present within a subset defined by
#' the provided arguments. It also verifies that there are no duplicate elements in the subset.
#' The function returns `TRUE` if all elements in the `whole_set` are found in the subset and
#' the subset contains only unique elements; otherwise, it returns `FALSE`.
#'
#' @param whole_set A vector containing the elements to be checked for membership in the subset.
#' @param ... Elements defining the subset. All provided elements are combined into a single vector
#'   and checked for uniqueness and membership.
#'
#' @return A logical value indicating whether all elements of `whole_set` are present in the subset
#'   and if the subset contains only unique elements. Returns `TRUE` if the conditions are met;
#'   otherwise, returns `FALSE`.
#'
#' @export
#'
#' @examples
#' # Example with all elements of whole_set in the subset
#' whole_set <- c(1, 2, 3)
#' subset <- c(1, 2, 3)
#' is.subset_membership(whole_set, subset) # Should return TRUE
#'
#' # Example with missing elements in the subset
#' whole_set <- c(1, 2, 3)
#' subset_1 <- c(1, 2)
#' subset_2 <- c(3)
#' is.subset_membership(whole_set, subset_1, subset_2) # Should return TRUE
#'
#' # Example with duplicate elements in the subset
#' whole_set <- c(1, 2, 3)
#' subset <- c(1, 2, 2, 3)
#' is.subset_membership(whole_set, subset) # Should return FALSE
is.subset_membership <- function(whole_set, ...) {
  subset <- c(...)

  # Check if all elements in whole_set are present in the subset and that subset contains only unique elements
  if (!all(whole_set %in% subset) || length(unique(subset)) != length(subset)) {
    return(FALSE)
  }

  return(TRUE)
}

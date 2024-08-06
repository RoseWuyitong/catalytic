# ---- GLM Functions ----

#' Calculate GLM Lambda Based on Discrepancy Method
#'
#' This function calculates the GLM lambda based on the specified discrepancy method.
#'
#' @param discrepancy_method Character. The method used to calculate the discrepancy.
#'   One of "square_error", "classification_error", or "logistic_deviance".
#' @param X A matrix of predictor variables.
#' @param coefs A numeric vector of coefficients.
#'
#' @return Numeric value of the calculated GLM lambda.
#'
#' @examples
#' X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
#' coefs <- c(0.5, 0.2)
#' get_glm_lambda("square_error", X, coefs)
#' get_glm_lambda("classification_error", X, coefs)
#' get_glm_lambda("logistic_deviance", X, coefs)
get_glm_lambda <- function(
    discrepancy_method = c(
      "square_error",
      "classification_error",
      "logistic_deviance"
    ),
    X,
    coefs) {
  discrepancy_method <- match.arg(discrepancy_method)

  switch(discrepancy_method,
    "square_error" = 2 * get_glm_mean("gaussian", X, coefs),
    "classification_error" = 2 * get_glm_mean("binomial", X, coefs),
    "logistic_deviance" = get_linear_predictor(X, coefs)
  )
}

#' Calculate Diagonal Approximate Covariance
#'
#' This function calculates the diagonal of the approximate covariance matrix for a GLM model.
#'
#' @param X A matrix of predictor variables.
#' @param model A GLM model object.
#'
#' @return A numeric vector containing the diagonal elements of the approximate covariance matrix.
#'
#' @examples
#' data(mtcars)
#' model <- glm(mpg ~ wt + hp, data = mtcars, family = binomial())
#' get_glm_diag_approx_cov(mtcars[, c("wt", "hp")], model)
get_glm_diag_approx_cov <- function(X, model) {
  # Function to compute the second derivative matrix
  get_second_derivative <- function(X,
                                    fitted_values,
                                    prior_weights) {
    return((t(cbind(1, X))) %*%
      diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
      cbind(1, X))
  }

  return(diag(chol2inv(chol(get_second_derivative(
    X = as.matrix(X),
    fitted_values = model$fitted.values,
    prior_weights = model$prior.weights
  )))))
}

#' Get GLM Family Name
#'
#' This function returns the name of a generalized linear model (GLM) family.
#' Optionally, it can also return the family name with its link function.
#'
#' @param family A character string naming a family function, a family function, or the result of a call to a family function.
#' @param with_link A logical value indicating whether to return the family name with its link function.
#'
#' @return A character string representing the family name, optionally including the link function.
#'
#' @examples
#' get_glm_family(gaussian)
#' get_glm_family("gaussian", with_link = TRUE)
#'
#' get_glm_family("binomial")
#' get_glm_family(binomial, with_link = TRUE)
get_glm_family <- function(family, with_link = FALSE) {
  validate_family_glm <- function(f) {
    if (is.character(f)) {
      tryCatch(
        {
          f <- get(f, mode = "function", envir = parent.frame(2))
        },
        error = function(.) stop("'family' must be a family.", call. = FALSE)
      )
    }
    if (is.function(f)) {
      f <- f()
    }
    if (!is(f, "family")) {
      stop("'family' must be a family.", call. = FALSE)
    }

    return(f)
  }

  family <- validate_family_glm(family)

  if (with_link) {
    return(paste0(family$family, " [", family$link, "]"))
  }

  return(family$family)
}

#' Calculate the Mean for a GLM
#'
#' This function calculates the mean of the response variable for a given generalized linear model (GLM) based on the family type, predictor variables, and coefficients.
#'
#' @param family A character string naming a family function, a family function, or the result of a call to a family function.
#' @param X A numeric matrix or data frame of predictor variables.
#' @param coefs A numeric vector of coefficients.
#'
#' @return A numeric vector representing the mean of the response variable.
#'
#' @examples
#' X <- matrix(1:6, ncol = 2)
#' coefs <- c(0.5, 0.5)
#' get_glm_mean(gaussian(), X, coefs)
#'
#' coefs_with_intercept <- c(1, 0.5, 0.5)
#' get_glm_mean("binomial", X, coefs_with_intercept)
get_glm_mean <- function(family,
                         X,
                         coefs) {
  lp <- get_linear_predictor(X, coefs)

  switch(get_glm_family(family),
    "binomial" = 1 / (1 + exp(-lp)),
    "gaussian" = lp,
    stop("Invalid GLM family. Supported families are 'binomial' and 'gaussian'.", call. = FALSE)
  )
}

#' Calculate the Log Density for a GLM
#'
#' This function calculates the log density of the response variable for a given generalized linear model (GLM) based on the family type, predictor variables, coefficients, and weights.
#'
#' @param family A character string naming a family function, a family function, or the result of a call to a family function.
#' @param X A numeric matrix or data frame of predictor variables.
#' @param Y A numeric vector of response variables.
#' @param coefs A numeric vector of coefficients.
#' @param weights A numeric vector of weights (optional, default is 1).
#'
#' @return A numeric value representing the log density of the response variable.
#'
#' @examples
#' X <- matrix(1:6, ncol = 2)
#' Y <- c(0, 1, 0)
#' coefs <- c(0.5, 0.5)
#' get_glm_log_density(binomial(), X, Y, coefs)
#'
#' coefs_with_intercept <- c(1, 0.5, 0.5)
#' get_glm_log_density("gaussian", X, Y, coefs_with_intercept)
get_glm_log_density <- function(family,
                                X,
                                Y,
                                coefs,
                                weights = 1) {
  lp <- get_linear_predictor(X, coefs)
  if (is.data.frame(Y)) {
    as.vector(apply(Y, 2, function(x) as.numeric(as.character(x))))
  } else {
    as.vector(Y)
  }

  switch(get_glm_family(family),
    "binomial" = sum(weights * (lp * Y - log(1 + exp(lp)))),
    "gaussian" = 0.5 * sum(weights * (1 + log(2 * pi * var(Y - lp)))),
    stop("Invalid GLM family. Supported families are 'binomial' and 'gaussian'.", call. = FALSE)
  )
}

#' Calculate the Gradient of Log Density for a GLM
#'
#' This function calculates the gradient of the log density of the response variable
#' for a given generalized linear model (GLM) based on the family type, predictor variables, coefficients, and weights.
#'
#' @param family A character string naming a family function, a family function, or
#'               the result of a call to a family function.
#' @param X A numeric matrix or data frame of predictor variables.
#' @param Y A numeric vector of response variables.
#' @param coefs A numeric vector of coefficients.
#' @param weights A numeric vector of weights (optional, default is 1).
#'
#' @return A numeric vector representing the gradient of the log density of the response variable.
#'
#' @examples
#' X <- matrix(1:6, ncol = 2)
#' Y <- c(0, 1, 0)
#' coefs <- c(0.5, 0.5)
#' get_log_density_grad_glm(binomial(), X, Y, coefs)
#'
#' coefs_with_intercept <- c(1, 0.5, 0.5)
#' get_log_density_grad_glm(gaussian(), X, Y, coefs_with_intercept)
get_glm_log_density_grad <- function(family,
                                     X,
                                     Y,
                                     coefs,
                                     weights = 1) {
  return(c(t(weights * (as.vector(Y) - get_glm_mean(family, X, coefs))) %*% as.matrix(cbind(1, X))))
}

#' Generate Simulated Data
#'
#' This function generates simulated data based on the specified family distribution.
#'
#' @param family The family distribution of the data. Defaults to 'gaussian'.
#' @param n The number of observations to generate. Defaults to 10.
#' @param mean The mean value for the distribution. Defaults to 0.
#' @param sd The standard deviation for the distribution. Defaults to 1.
#'
#' @return A numeric vector of simulated data.
#'
#' @examples
#' # Generate 10 observations from a Gaussian distribution
#' get_sample_data()
#'
#' # Generate 20 observations from a binomial distribution with mean 0.5
#' get_sample_data(family = "binomial", n = 20, mean = 0.5)
get_glm_sample_data <- function(family = gaussian,
                                n = 10,
                                mean = 0,
                                sd = 1) {
  switch(get_glm_family(family),
    "binomial" = stats::rbinom(
      n,
      1,
      mean
    ),
    "gaussian" = stats::rnorm(
      n,
      mean,
      sd
    )
  )
}


give_glm_bayes_joint_suggestion <- function(family,
                                            stan_iter,
                                            alpha,
                                            stan_sample_model,
                                            binomial_use_theta,
                                            binomial_joint_alpha,
                                            binomial_tau_lower) {
  get_est_params <- function(params, percentiles, mean_, distribution) {
    alpha <- exp(params[1])
    beta <- exp(params[2])
    calc_percentiles <- switch(distribution,
      "invgamma" = invgamma::qinvgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta),
      "gamma" = stats::qgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta)
    )
    error <- sum((percentiles - calc_percentiles)^2) + (mean_ - (beta / (alpha - 1)))^2
    return(error)
  }

  if (get_glm_family(family) == "binomial") {
    if (binomial_use_theta) {
      theta_summary <- rstan::summary(stan_sample_model)$summary["theta", ]

      est_params <- stats::optim(
        par = c(1, 1),
        fn = get_est_params,
        percentiles = theta_summary[c("25%", "50%", "75%")],
        mean_ = theta_summary["mean"],
        distribution = "invgamma"
      )

      if (invgamma::qinvgamma(c(0.975),
        shape = exp(est_params$par[1]),
        scale = exp(est_params$par[2])
      ) < theta_summary["97.5%"]) {
        if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
          cat("\nSUGGESTION: Please consider using more iteration for stan (e.g.", stan_iter + 1000, ")\n")
        } else if (alpha < 5) {
          cat("\nSUGGESTION: Please consider using larger `alpha` (e.g.", min(alpha + 0.1, 5), ")\n")
        } else {
          cat("\nSUGGESTION: Please consider add `binomial_joint_alpha = TRUE`.\n")
        }
      }
    } else {
      tau_summary <- rstan::summary(stan_sample_model)$summary["tau", ]

      est_params <- stats::optim(
        par = c(1, 1),
        fn = get_est_params,
        percentiles = tau_summary[c("25%", "50%", "75%")],
        mean_ = tau_summary["mean"],
        distribution = "gamma"
      )


      if (stats::qgamma(c(0.025),
        shape = exp(est_params$par[1]),
        scale = exp(est_params$par[2])
      ) > tau_summary["2.5%"]) {
        if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
          cat("\nSUGGESTION: Please consider using more iteration for stan (e.g.", stan_iter + 1000, ")\n")
        } else if (binomial_tau_lower < 1) {
          cat("\nSUGGESTION: Please consider using larger `binomial_tau_lower` (e.g.", min(binomial_tau_lower + 0.05, 1), ")\n")
        } else {
          cat("\nSUGGESTION: Please consider using `binomial_use_theta = TRUE`.\n")
        }
      }
    }
  }
}

#' Check Dataset Size and Risk Estimate Method
#'
#' This function checks the size of the observational dataset in relation to the specified risk estimate method.
#' It provides warnings and prompts based on the dataset size to guide the user towards optimal methods.
#'
#' @param obs_size An integer representing the size of the observational dataset.
#' @param risk_estimate_method A character string specifying the risk estimate method.
#'                             Possible values are "parametric_bootstrap", "steinian_estimate", and "cross_validation".
#'
#' @return None
check_glm_tune_size_and_risk_estimate_method <- function(obs_size,
                                                         risk_estimate_method) {
  if (obs_size > if (risk_estimate_method == "parametric_bootstrap") {
    600
  } else if (risk_estimate_method == "steinian_estimate") {
    400
  } else {
    Inf
  }) {
    response <- menu(c("YES", "NO"), title = paste0(
      "The large size of the dataset ",
      "could make estimate time longer than expected.",
      "\nDo you want to continue?\n"
    ))

    if (response == 2) {
      stop("Process terminated by user.", call. = FALSE)
    }
  }


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

  if (obs_size < if (risk_estimate_method == "cross_validation") 30 else 0) {
    warning(paste0(
      "The size of the observational dataset is insufficient ",
      "for cross vlidation to yield optimal results.",
      "\nSUGGESTION: Please consider utilizing alternative methods",
      "(e.g. risk_estimate_method = parametric_bootstrap) ",
      "for better estimate."
    ), call. = FALSE)
  }
}

# ---- COX Functions ----


#' Solve a linear system using QR decomposition
#'
#' This function solves a linear system of equations represented by `coef_matrix %*% x = constant_vector`
#' using QR decomposition. It handles NA values in the solution by replacing them with a small constant.
#'
#' @param coef_matrix A matrix of coefficients from the linear system.
#' @param constant_vector The constant vector from the linear system.
#'
#' @return The solution vector `x` from the linear system.
#'
#' @examples
#' coef_matrix <- matrix(c(1, 1, 1, 2, 3, 5, 4, 6, 8), nrow = 3, byrow = TRUE)
#' constant_vector <- c(10, 20, 30)
#' solution <- get_cox_qr_solve(coef_matrix, constant_vector)
#' print(solution)
get_cox_qr_solve <- function(coef_matrix,
                             constant_vector) {
  # Calculate solution using QR decomposition
  qr_solution <- solve(coef_matrix, constant_vector)

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
}

#' Compute the Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function calculates the Hessian matrix for a Cox proportional hazards model,
#' which is used to estimate the second-order partial derivatives of the log-likelihood function.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-status or censoring time for each observation.
#' @param status A numeric vector indicating whether the status of interest occurred (1) or the observation was censored (0).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector of entry times for each observation. Default is NULL, in which case all entry times are set to 0.
#'
#' @return A numeric matrix representing the Hessian matrix.
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' status <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
#' coefs <- rnorm(ncol(X))
#'
#' # Compute Hessian matrix
#' hessian_matrix <- get_cox_hessian(X, time, status, coefs)
#' print(hessian_matrix)
get_cox_hessian <- function(X,
                            time,
                            status,
                            coefs,
                            entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  cersored_status <- which(as.numeric(status) == 1)
  hessian <- matrix(0, ncol(X), ncol(X))

  for (i in cersored_status) {
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    risk_set <- matrix(as.matrix(X[risk_set_idx, ]), nrow = length(risk_set_idx))
    exp_risk_lp <- exp(risk_set %*% coefs)
    sum_exp_risk_lp <- sum(exp_risk_lp)
    risk_exp_risk_lp <- t(risk_set) %*% exp_risk_lp

    risk_sqrt_exp_risk_lp <- sapply(
      1:length(risk_set_idx),
      function(i) (t(risk_set)[, i] * sqrt(exp_risk_lp[i]))
    )

    hessian <- hessian +
      tcrossprod(as.matrix(risk_sqrt_exp_risk_lp)) / sum_exp_risk_lp -
      tcrossprod(as.matrix(risk_exp_risk_lp)) / sum_exp_risk_lp^2
  }
  return(-hessian)
}

#' Compute the Synthetic Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function calculates the synthetic Hessian matrix for a Cox proportional hazards model,
#' using the provided predictor variables, time-to-status data, coefficients, and a hazard constant.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-status or censoring time for each observation.
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param hazard_constant A numeric constant representing the hazard.
#'
#' @return A numeric matrix representing the synthetic Hessian matrix.
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' coefs <- rnorm(ncol(X))
#' hazard_constant <- 0.1
#'
#' # Compute synthetic Hessian matrix
#' hessian_matrix <- get_cox_syn_gradient(X, time, coefs, hazard_constant)
#' print(hessian_matrix)
get_cox_syn_hessian <- function(X,
                                time,
                                coefs,
                                hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)

  syn_lp <- as.numeric(get_linear_predictor(X, coefs))

  hessian <- matrix(0, ncol(X), ncol(X))

  for (i in 1:nrow(X)) {
    hessian <- hessian -
      hazard_constant *
        time[i] *
        exp(syn_lp[i]) *
        tcrossprod(X[i, ])
  }

  return(hessian)
}

#' Compute the Gradient for Cox Proportional Hazards Model
#'
#' This function calculates the gradient for a Cox proportional hazards model using the provided predictor variables,
#' time-to-status data, status indicators, coefficients, and optional entry points.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-status or censoring time for each observation.
#' @param status A numeric vector of status indicators (1 if status occurred, 0 if censored).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector of entry points for each observation. Default is NULL.
#'
#' @return A numeric vector representing the gradient for the Cox proportional hazards model.
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' status <- c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1)
#' coefs <- rnorm(ncol(X))
#'
#' # Compute gradient
#' gradient <- get_cox_gradient(X, time, status, coefs)
#' print(gradient)
get_cox_gradient <- function(X,
                             time,
                             status,
                             coefs,
                             entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  gradient <- apply(
    cbind(time, X), 1,
    function(df) {
      df <- matrix(df, nrow = 1)
      ts <- df[, 1]
      xs <- df[, 2:ncol(df)]
      X_risk_set <- as.matrix(X[get_cox_risk_set_idx(ts, entry_points, time, status), ])

      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% coefs)
      t2 <- sum(exp(X_risk_set %*% coefs))

      return(xs - t1 / t2)
    }
  ) %*% status

  return(as.vector(gradient))
}

#' Compute the Gradient for Synthetic Data in a Cox Proportional Hazards Model
#'
#' This function calculates the gradient for synthetic data in a Cox proportional hazards model
#' using the provided predictor variables, time-to-status data, coefficients, and a hazard constant.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-status or censoring time for each observation.
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param hazard_constant A numeric value representing the hazard constant.
#'
#' @return A numeric vector representing the gradient for the synthetic data in the Cox proportional hazards model.
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' coefs <- rnorm(ncol(X))
#' hazard_constant <- 0.01
#'
#' # Compute gradient
#' gradient <- get_cox_syn_gradient(X, time, coefs, hazard_constant)
#' print(gradient)
get_cox_syn_gradient <- function(X,
                                 time,
                                 coefs,
                                 hazard_constant) { # TODO Optimize code, e.g. replace for loop

  X <- as.matrix(X)
  time <- as.vector(time)

  syn_lp <- as.numeric(get_linear_predictor(X, coefs))
  gradient <- rep(0, ncol(X))
  for (i in 1:nrow(X)) {
    gradient <- gradient + (1 - hazard_constant * time[i] * exp(syn_lp[i])) * X[i, ]
  }

  return(as.numeric(gradient))
}


#' Compute the Partial Likelihood for a Cox Proportional Hazards Model
#'
#' This function calculates the partial likelihood for a Cox proportional hazards model
#' using the provided predictor variables, time-to-status data, status indicators, and coefficients.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-status or censoring time for each observation.
#' @param status A numeric vector indicating whether an status occurred (1) or the data is censored (0).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector representing the entry times for each observation. Default is NULL.
#'
#' @return A numeric value representing the partial likelihood for the Cox proportional hazards model.
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' status <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 1)
#' coefs <- rnorm(ncol(X))
#' entry_points <- NULL
#'
#' # Compute partial likelihood
#' pl <- get_cox_partial_likelihood(X, time, status, coefs, entry_points)
#' print(pl)
get_cox_partial_likelihood <- function(X,
                                       time,
                                       status,
                                       coefs,
                                       entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  cersored_status <- which(as.numeric(status) == 1)
  pl <- 0
  for (i in cersored_status) {
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    exp_risk_lp <- exp(as.vector(get_linear_predictor(
      matrix(as.matrix(X[risk_set_idx, ]), nrow = length(risk_set_idx)),
      coefs
    )))
    pl <- pl + sum(X[i, ] * coefs) - log(sum(exp_risk_lp))
  }

  return(pl)
}


get_cox_kappa <- function(X,
                          time,
                          status,
                          hazard_constant = NULL,
                          tol = 1e-5,
                          max_iter = 25,
                          entry_points = NULL,
                          initial_values = NULL) {
  if (is.null(hazard_constant)) hazard_constant <- sum(status) / mean(time) / nrow(X)
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  coefs <- if (is.null(initial_values)) rep(0, ncol(X)) else initial_values

  iter <- 1
  diff_new_old <- TRUE
  while ((iter < max_iter) & diff_new_old) {
    coefs_old <- coefs

    exp_lp <- exp(get_linear_predictor(X, coefs_old))
    gd_old <- get_cox_syn_gradient(
      X = X,
      time = time,
      coefs = coefs_old,
      hazard_constant = hazard_constant
    )
    hess_old <- get_cox_syn_hessian(
      X = X,
      time = time,
      coefs = coefs_old,
      hazard_constant = hazard_constant
    )

    coefs <- coefs_old - get_cox_qr_solve(hess_old, gd_old)

    iter <- iter + 1

    if (any(is.na(coefs))) stop("Please note NA appears in iteration.")

    if (norm(as.numeric(coefs), "2") > 1000000) stop("Please not the euclidean norm of coefficients is too large.")

    diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > tol * norm(as.numeric(coefs_old), "2")
  }

  lp_new <- as.numeric(get_linear_predictor(X, coefs))

  return(sum(lp_new - exp(lp_new) * hazard_constant * time) / nrow(X))
}

#' Compute Risk and Failure Sets for Survival Analysis
#'
#' This function computes the risk and failure sets for a given time vector, status vector, and a sequence of times (indicator_vector).
#' It identifies the times at which statuss occur or subjects are censored.
#'
#' @param time_vector A numeric vector representing the time-to-status or censoring time for each observation.
#' @param status_vector A numeric vector indicating whether an status occurred (1) or the data is censored (0).
#' @param indicator_vector A numeric vector of times at which to evaluate the risk and failure sets.
#'
#' @return A list containing two matrices: `risk_set` and `failure_set`, representing the risk and failure sets at each time point.
#'
#' @examples
#' time_vector <- c(5, 8, 12, 18, 22)
#' status_vector <- c(1, 0, 1, 1, 0)
#' indicator_vector <- seq(1, 25, by = 1)
#' result <- get_risk_and_failure_sets(time_vector, status_vector, s)
#' print(result$risk_set)
#' print(result$failure_set)
get_cox_risk_and_failure_sets <- function(time_vector,
                                          status_vector,
                                          indicator_vector) {
  time_vector <- as.vector(time_vector)
  status_vector <- as.vector(status_vector)
  indicator_vector <- as.vector(indicator_vector)

  # Function to find the time where the status occurred or reach censoring time
  find_occur_idx <- function(time) {
    min(which(indicator_vector >= time))
  }

  risk_set <- failure_set <- matrix(
    0,
    length(time_vector),
    length(indicator_vector)
  )

  # Process subjects with status
  for (i in which(status_vector == 1)) {
    # Mark only exactly when the status occurrence
    failure_set[i, find_occur_idx(time_vector[i])] <- 1
    # Mark from risk up to the status time
    risk_set[i, 1:find_occur_idx(time_vector[i])] <- 1
  }

  # Process subjects without status
  for (i in which(status_vector == 0)) {
    # Mark the subject as at risk up to the censoring time
    risk_set[i, 1:find_occur_idx(time_vector[i])] <- 1
  }

  return(list(
    risk_set = risk_set,
    failure_set = failure_set
  ))
}

#' Get Risk Set Index for Survival Analysis
#'
#' This function identifies the risk set for a given time of interest in survival analysis.
#' It returns the indices of subjects who are at risk at the specified time.
#'
#' @param time_of_interest A numeric value representing the time of interest.
#' @param entry_vector A numeric vector representing the entry times for each subject.
#' @param time_vector A numeric vector representing the time-to-status or censoring time for each subject.
#' @param status_vector A numeric vector indicating whether an status occurred (1) or the data is censored (0).
#'
#' @return A numeric vector of indices representing the subjects at risk at the specified time of interest.
#'
#' @examples
#' time_of_interest <- 10
#' entry_vector <- c(0, 0, 0, 5, 7)
#' time_vector <- c(5, 8, 12, 18, 22)
#' status_vector <- c(1, 0, 1, 1, 0)
#' risk_set_idx <- get_cox_risk_set_idx(time_of_interest, entry_vector, time_vector, status_vector)
#' print(risk_set_idx)
get_cox_risk_set_idx <- function(time_of_interest,
                                 entry_vector,
                                 time_vector,
                                 status_vector) {
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & status_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}

# ---- LMM Functions ----

#' Update Sigma Square for LMM
#'
#' This function performs the minimization program for updating the sigma square
#' parameter in a linear mixed model (LMM). It computes the sum of the
#' log of the variances plus the squared terms divided by the variances, adjusted by tau.
#'
#' @param sigma_square Numeric value representing the sigma square parameter.
#' @param sigma_u_square Numeric value representing the sigma_u square parameter.
#' @param d_seq Numeric vector of d sequence values.
#' @param d_seq_star Numeric vector of d star sequence values.
#' @param q Numeric vector of q values.
#' @param q_star Numeric vector of q star values.
#' @param tau Numeric value for the tau parameter.
#'
#' @return Numeric value representing the result of the minimization program.
#'
#' @examples
#' # Example data
#' sigma_square <- 0.5
#' sigma_u_square <- 1.0
#' d_seq <- c(1, 2, 3)
#' d_seq_star <- c(1.5, 2.5, 3.5)
#' q <- c(0.1, 0.2, 0.3)
#' q_star <- c(0.15, 0.25, 0.35)
#' tau <- 0.01
#'
#' # Update sigma square
#' result <- update_sigma_square_cat(sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau)
#' print(result)
update_lmm_sigma_square_cat <- function(sigma_square,
                                        sigma_u_square,
                                        d_seq,
                                        d_seq_star,
                                        q,
                                        q_star,
                                        tau) {
  return(sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) +
    tau / length(d_seq_star) * sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square)))
}


# ---- Other Functions ----

get_h_matrix <- function(obs_x, syn_x, tau) {
  obs_x <- as.matrix(obs_x)
  syn_x <- as.matrix(syn_x)
  obs_size <- nrow(obs_x)
  syn_size <- nrow(syn_x)

  return(obs_x %*% (
    chol2inv(chol((t(obs_x) %*% obs_x + (tau / syn_size) * (t(syn_x) %*% syn_x))))) %*% (
    t(obs_x) + tau / (syn_size * obs_size) * t(syn_x) %*% rep(1, syn_size) %*% t(rep(1, obs_size))
  ))
}


#' Memoized Function to Retrieve Stan Models
#'
#' This function memoizes the creation of a Stan model based on the specified model type
#' and optional parameters for specific methods like GLMs. It returns a compiled Stan model object.
#'
#' @param type A character string specifying the type of Stan model, choose from c("glm", "cox").
#' @param glm_family A character string specifying the family for GLM models. Default is NULL.
#' @param joint_tau Logical indicating joint modeling for tau parameter. Default is FALSE.
#' @param glm_binomial_joint_theta Logical indicating joint modeling for theta parameter in binomial GLMs. Default is FALSE.
#' @param glm_binomial_joint_alpha Logical indicating joint modeling for alpha parameter in binomial GLMs. Default is FALSE.
#' @param glm_binomial_gibbs Logical indicating whether Gibbs sampling is used in binomial GLMs. Default is FALSE.
#'
#' @return A compiled Stan model object.
#'
#' @examples
#' # Retrieve a Stan model for binomial logistic regression with joint theta
#' model <- get_stan_model("glm", glm_family = "binomial", glm_binomial_joint_theta = TRUE)
get_stan_model <- memoise::memoise(function(type = c("glm", "cox"),
                                            glm_family = NULL,
                                            joint_tau = FALSE,
                                            glm_binomial_joint_theta = FALSE,
                                            glm_binomial_joint_alpha = FALSE,
                                            glm_binomial_gibbs = FALSE) {
  glm_family <- if (!is.null(glm_family)) get_glm_family(glm_family) else NULL

  return(rstan::stan_model(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      paste0(
        c(
          type,
          glm_family,
          if (any(
            joint_tau,
            glm_binomial_joint_theta,
            glm_binomial_joint_alpha
          )) {
            "joint"
          },
          if (type == "glm" && glm_family == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs) "theta",
          if (type == "glm" && glm_family == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs & glm_binomial_joint_alpha) "alpha"
        ),
        collapse = "_"
      ),
      ".stan"
    )
  ))
})


#' Adjust Data Based on Formula
#'
#' This function adjusts the given data according to the provided formula. It extracts the model frame
#' from the formula and applies it to the dataset, separating the adjusted data into observed and synthetic parts.
#'
#' @param cat_init A list from cat_glm_initialization, cat_lmm_initialization or cat_cox_initialization
#' @param formula A formula specifying the adjustment to be applied to the `x` matrix.
#'
#' @return The input list with additional elements:
#'   \item{`adj_x`} {The adjusted data matrix.}
#'   \item{`adj_obs_x`} {The observed part of the adjusted data.}
#'   \item{`adj_syn_x`} {The synthetic part of the adjusted data.}
get_adjusted_data <- function(cat_init, formula) {
  cat_init$adj_x <- model.frame(
    get_formula_rhs(as.formula(formula), TRUE),
    as.data.frame(cat_init$x)
  )

  cat_init$adj_x <- as.data.frame(apply(cat_init$adj_x, 2, function(x) as.numeric(as.character(x))))

  rownames(cat_init$adj_x) <- rownames(cat_init$x)

  cat_init$adj_obs_x <- as.data.frame(cat_init$adj_x[1:cat_init$obs_size, ])
  cat_init$adj_syn_x <- as.data.frame(cat_init$adj_x[(cat_init$obs_size + 1):cat_init$size, ])

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

#' Apply Resampling Methods to Data
#'
#' This function applies various resampling techniques to a given vector \code{x} based on the
#' original data \code{original_x}. The available methods include symmetrizing, smoothing,
#' flattening, and deskewing.
#'
#' @param x A numeric vector to be resampled.
#' @param original_x The original numeric vector from which \code{x} was derived.
#' @param method A character string specifying the resampling method to be used.
#'               Options are "symmetrizing", "smoothing", "flattening", and "deskewing".
#'               Default is "symmetrizing".
#'
#' @return A numeric vector with resampled values based on the specified method.
#'
#' @details
#' \strong{symmetrizing}: Adjusts \code{x} by adding random variations around the mean of \code{original_x}.
#'
#' \strong{smoothing}: Adds random noise to \code{x} based on the range of values in \code{original_x}.
#'
#' \strong{flattening}: Randomly replaces half of the values in \code{x} with values sampled from \code{original_x}.
#'
#' \strong{deskewing}: Adjusts the skewness of \code{x} by sampling values from a truncated normal distribution
#'                   based on \code{original_x}.
#'
#' @examples
#' x <- rnorm(10)
#' original_x <- rnorm(10)
#'
#' # Apply symmetrizing method
#' resampled_x_sym <- get_resampling_methods(x, original_x, method = "symmetrizing")
#' print(resampled_x_sym)
#'
#' # Apply smoothing method
#' resampled_x_smooth <- get_resampling_methods(x, original_x, method = "smoothing")
#' print(resampled_x_smooth)
#'
#' # Apply flattening method
#' resampled_x_flat <- get_resampling_methods(x, original_x, method = "flattening")
#' print(resampled_x_flat)
#'
#' # Apply deskewing method
#' resampled_x_deskew <- get_resampling_methods(x, original_x, method = "deskewing")
#' print(resampled_x_deskew)
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

  method <- match.arg(method)

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

#' Standardize and preprocess a data frame
#'
#' This function converts columns in a data frame to numeric or factor types,
#' and replaces missing values with either a specified value or the mode of the column.
#'
#' @param data A data frame containing numeric or categorical data.
#' @param na_replace Function or value to replace NA values. Default is `stats::na.omit`.
#'
#' @return A standardized data frame with converted columns and replaced NA values.
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, NA, 4),
#'   x2 = c("A", "B", NA, "A")
#' )
#' standardized_data_omit <- get_standardized_data(data)
#' print(standardized_data_omit)
#'
#' standardized_data_mean <- get_standardized_data(data, mean)
#' print(standardized_data_mean)
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

#' Calculate the Linear Predictor
#'
#' This function calculates the linear predictor (linear combination of predictors
#' and coefficients) for a given set of predictors and coefficients.
#'
#' @param X A numeric matrix or data frame of predictor variables.
#' @param coefs A numeric vector of coefficients.
#'              If the length of \code{coefs} is one more than the number of columns in \code{X},
#'              the first element of \code{coefs} is treated as the intercept.
#'
#' @return A numeric vector representing the linear predictor.
#'
#' @examples
#' X <- matrix(1:6, ncol = 2)
#' coefs <- c(0.5, 0.5)
#' get_linear_predictor(X, coefs)
#'
#' coefs_with_intercept <- c(1, 0.5, 0.5)
#' get_linear_predictor(X, coefs_with_intercept)
get_linear_predictor <- function(X, coefs) {
  # X <- apply(X, 2, function(x) as.numeric(as.character(x)))
  X <- as.matrix(X)
  coefs <- as.numeric(coefs)
  if (length(coefs) == ncol(X)) {
    return(as.vector(X %*% coefs))
  } else if (length(coefs) == (ncol(X) + 1)) {
    return(as.vector(coefs[1] + X %*% coefs[-1]))
  }

  stop(
    "Dimensions of coefs and predictor variables do not match.\n\n Input Coefficients:\n",
    coefs,
    "\n\nInput Covariance:\n",
    X
  )
}

#' Get Formula String Representation
#'
#' Converts a formula object into a compact string representation, removing excess whitespace.
#'
#' @param formula A formula object.
#'
#' @return A character string representing the formula.
#'
#' @examples
#' f <- y ~ x1 + x2
#' get_formula_string(f)
#'
#' g <- z ~ w
#' get_formula_string(g)
get_formula_string <- function(formula) {
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  return(char)
}

#' Get Formula Left-hand Side (LHS) Representation
#'
#' Extracts the left-hand side of a formula as a compact string representation.
#'
#' @param formula A formula object.
#'
#' @return A character string representing the left-hand side of the formula.
#'
#' @examples
#' f <- y ~ x1 + x2
#' get_formula_lhs(f)
#'
#' g <- z ~ w
#' get_formula_lhs(g)
get_formula_lhs <- function(formula) {
  return(get_formula_string(rlang::f_lhs(formula)))
}


#' Get Formula Right-hand Side (RHS) Representation
#'
#' Extracts the right-hand side of a formula as a compact string representation.
#'
#' @param formula A formula object.
#' @param with_tilde Logical; if TRUE, includes the tilde (~) symbol before the RHS (default is FALSE).
#'
#' @return A character string representing the right-hand side of the formula.
#'
#' @examples
#' f <- y ~ x1 + x2
#' get_formula_rhs(f)
#'
#' g <- z ~ w
#' get_formula_rhs(g, with_tilde = TRUE)
get_formula_rhs <- function(formula, with_tilde = FALSE) {
  if (with_tilde) {
    return(paste("~", deparse(rlang::f_rhs(formula)), collapse = " "))
  }
  return(get_formula_string(rlang::f_rhs(formula)))
}

#' Compute Discrepancy Measures for GLM Models
#'
#' Computes various discrepancy measures between observed and estimated responses
#' based on the specified discrepancy method.
#'
#' @param discrepancy_method Character string specifying the discrepancy method to use.
#'   Supported methods are "logarithmic_error", "square_error", "classification_error", and "logistic_deviance".
#' @param family GLM family specification (e.g., "binomial" or "gaussian").
#' @param X Predictor variables (matrix or data frame).
#' @param Y Observed response variable (numeric vector for square_error and logistic_deviance,
#'   factor or binary for classification_error).
#' @param coefs Coefficients vector from the GLM model.
#' @param est_Y Estimated response variable
#'
#' @return Numeric value representing the computed discrepancy measure.
#'
#' @examples
#' family <- "binomial"
#' X <- matrix(rnorm(10), ncol = 2)
#' Y <- rbinom(5, 1, 0.5)
#' coefs <- c(-1, 0.5, 0.3) # Example coefficients
#'
#' get_discrepancy("classification_error", family, X, Y, coefs)
#' get_discrepancy("logistic_deviance", family, X, Y, coefs)
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
  discrepancy_method <- match.arg(discrepancy_method)
  Y <- if (is.data.frame(Y)) {
    as.vector(apply(Y, 2, function(x) as.numeric(as.character(x))))
  } else {
    as.vector(Y)
  }

  if (is.null(est_Y) & !is.null(family)) est_Y <- get_glm_mean(family, X, coefs)

  if (discrepancy_method == "logarithmic_error") {
    Y <- pmax(0.0001, pmin(0.9999, Y))
    est_Y <- pmax(0.0001, pmin(0.9999, est_Y))
  }

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

#' Check if a List is Continuous
#'
#' This function checks if a given list is continuous. A list is considered
#' continuous if it is numeric and contains more than two unique values.
#'
#' @param lst A list or vector to be checked.
#'
#' @return TRUE if the list is continuous, FALSE otherwise.
#'
#' @examples
#' # Should return TRUE
#' is.continuous(c(1, 2, 3))
#'
#' # Should return FALSE
#' is.continuous(c(1, 1, 1))
#' is.continuous(c("a", "b", "c"))
is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}


#' Check if an object is a formula.
#'
#' This function checks whether the provided object \code{f} is a formula by
#' searching for the tilde (\code{~}) symbol in its deparsed form.
#'
#' @param f An object to be checked.
#'
#' @return \code{TRUE} if \code{f} is a formula (contains '~'), otherwise \code{FALSE}.
#'
#' @examples
#' # should return TRUE
#' is.formula(~ x + y)
#' is.formula(y ~ x)
#'
#' # should return FALSE
#' is.formula("x + y")
is.formula <- function(f) {
  return(grepl("~", deparse(f)))
}

#' Check if an object is a GLM formula.
#'
#' This function checks whether the provided object \code{f} is a Generalized Linear Model (GLM) formula.
#' It first verifies if \code{f} is not a formula for Linear Mixed Models (LMM) or Survival Analysis Models,
#' in which case it returns \code{FALSE}. Otherwise, it checks if \code{f} is a formula using the
#' \code{is.formula} function from the stats package.
#'
#' @param f An object to be checked.
#'
#' @return \code{TRUE} if \code{f} is a GLM formula and not an LMM or Survival formula, otherwise \code{FALSE}.
#'
#' @examples
#' # should return TRUE
#' is.formula_glm(~ x + y)
#' is.formula_glm(y ~ x)
#'
#' # should return FALSE
#' is.formula_glm(y ~ x + (a | b))
is.formula_glm <- function(f) {
  if (is.formula_lmm(f) | is.formula_surv(f)) {
    return(FALSE)
  }
  return(is.formula(f))
}

#' Check if an object is a formula for Linear Mixed Models (LMM).
#'
#' This function checks whether the provided object \code{f} is a formula specifically
#' designed for Linear Mixed Models (LMM). It does so by searching for patterns in the
#' deparsed form of \code{f}, specifically looking for expressions containing parentheses
#' (indicative of random effects in LMM formulas).
#'
#' @param f An object to be checked.
#'
#' @return \code{TRUE} if \code{f} is a formula for Linear Mixed Models (LMM), otherwise \code{FALSE}.
#'
#' @examples
#' # should return TRUE
#' is.formula_lmm(~ 1 + (1 | subject))
#' is.formula_lmm(Reaction ~ Days + (Days | Subject))
#'
#' # should return FALSE
#' is.formula_lmm(~ x + y)
#' is.formula_lmm("x + y")
is.formula_lmm <- function(f) {
  return(grepl("\\(.*|.*\\)", deparse(f)))
}

#' Check if an object is a survival formula.
#'
#' This function checks whether the provided object \code{f} is a survival formula
#' by searching for the occurrence of "Surv" in its deparsed form.
#'
#' @param f An object to be checked.
#'
#' @return \code{TRUE} if \code{f} is a survival formula (contains 'Surv'), otherwise \code{FALSE}.
#'
#' @examples
#' # Should return TRUE
#' is.formula_surv(Surv(time, status) ~ x)
#'
#' # Should return FALSE
#' is.formula_surv(~ x + y)
#' is.formula_surv("x + y")
is.formula_surv <- function(f) {
  return(grepl("Surv", deparse(f)))
}


#' Checks if the data object is generated from 'data_initilization_glm'
#'
#' This function checks if the given data object is generated from' data_initilization_glm'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_glm', otherwise FALSE.
is.cat_glm_init <- function(x) {
  return(inherits(x, "cat_initialization") & x$function_name == "cat_glm_initialization")
}

#' Checks if the data object is generated from 'data_initilization_lmm'
#'
#' This function checks if the given data object is generated from' data_initilization_lmm'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_lmm', otherwise FALSE.
is.cat_lmm_init <- function(x) {
  return(inherits(x, "cat_initilization") & x$function_name == "cat_lmm_initialization")
}

#' Checks if the data object is generated from 'data_initilization_cox'
#'
#' This function checks if the given data object is generated from' data_initilization_cox'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_cox', otherwise FALSE.
is.cat_cox_init <- function(x) {
  return(inherits(x, "cat_initialization") & x$function_name == "cat_cox_initialization")
}

#' Check if a set is a subset membership of another set.
#'
#' This function checks if all elements in \code{whole_set} are present in the
#' subset formed by \code{...} without overlap.
#'
#' @param whole_set A vector containing elements to be checked.
#' @param ... Vectors containing subsets to be checked.
#'
#' @return \code{TRUE} if \code{whole_set} is a subset of the union of all subsets,
#'         otherwise \code{FALSE}.
#'
#' @examples
#' # Should return TRUE
#' is.subset_membership(c(1, 2, 3), c(1, 2), c(3))
#'
#' # Should return FALSE
#' is.subset_membership(c(1, 2, 3), c(1, 4))
#' is.subset_membership(c(1, 2, 3), c(1, 2), c(4, 5))
is.subset_membership <- function(whole_set, ...) {
  subset <- c(...)

  # Check if there is any overlap
  if (length(unique(subset)) != sum(lengths(list(...)))) {
    return(FALSE)
  }

  # Check if all elements in whole_set are present in the subset
  if (!all(whole_set %in% subset)) {
    return(FALSE)
  }

  return(TRUE)
}

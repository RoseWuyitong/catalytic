# ----- Functions for data_initialization -----

#' Standardize and Clean Data Frame
#'
#' This function standardizes a data frame by replacing `NA` values and converting columns to appropriate data types.
#' It ensures that numeric columns remain numeric, logical columns remain logical, and other columns are converted to factors.
#'
#' @param data A data frame containing the data to be standardized and cleaned.
#' @param na_replace A value or function used to replace `NA` values.
#'                   If a numeric value is provided, it will replace all `NA` values.
#'                   If a function is provided (e.g., `mean` or `median`), the function will
#' .                  be applied to each column to replace `NA` values. Default is `stats::na.omit`.
#'
#' @return A standardized and cleaned data frame.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   continue_col = c(1, 2, NA, 4),
#'   bool_col =     c(1, 0, NA, 1),
#'   factor_col =   c("a", "b", "c", NA)
#' )
#'
#' # Example 1: Remove rows contains NA
#' standardized_data <- get_standardized_data(data)
#' print(standardized_data)
#'
#' # Example 2: Replace NA with a specific value and mode for factor columns
#' standardized_data <- get_standardized_data(data, na_replace = 0)
#' print(standardized_data)
#'
#' # Example 3: Replace NA with mean for numeric columns and mode for factor columns
#' standardized_data <- get_standardized_data(data, na_replace = base::mean)
#' print(standardized_data)
get_standardized_data <- function(data, na_replace = stats::na.omit) {
  convert_to_numeric_or_factor <- function(x) {
    if (is.continuous(stats::na.omit(x))) {
      x <- as.numeric(x)
    } else {
      x <- as.factor(x)
    }
  }

  replace_df_na <- function(x) {
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

  data <- as.data.frame(lapply(data, convert_to_numeric_or_factor))
  data <- as.data.frame(if (identical(stats::na.omit, na_replace)) {
    stats::na.omit(data)
  } else {
    lapply(data, replace_df_na)
  })

  if (nrow(data) == 1) {{ stop(
    "The data size is 1, which is too small. cannot process."
  ) }}

  return(data)
}

#' Resample Data Frame Columns
#'
#' This function performs resampling on each column of a data frame based on its type (continuous or categorical).
#' It provides options for various resampling methods such as smoothing for continuous variables and flattening for categorical variables.
#'
#' @param data A data frame containing the data to be resampled.
#' @param resample_size The size of the resampled data.
#' @param data_degree Optional vector specifying the degrees of each column in `data`.
#'                    If not provided, defaults to all columns having degree 1.
#' @param resample_only Logical indicating whether to return only the resampled data without recording the process. Default is FALSE.
#'
#' @return A list containing:
#'   - \code{resampled_data}: The resampled data frame.
#'   - \code{resampled_data_record}: A data frame recording the resampling process for each column.
#' @export
#'
#' @examples
#' # Example 1: Resample data frame columns with default settings
#' data <- data.frame(
#'   col1 = c(1, 1, 1, 1, 10),
#'   col2 = c("A", "B", "C", "D", "E")
#' )
#' resampled <- get_resampled_data(data, resample_size = 10)
#' print(resampled$resampled_data)
#' print(resampled$resampled_data_record)
#'
#' # Example 2: Resample data frame columns with specified resample size and data degree
#' resampled <- get_resampled_data(data, resample_size = 10, data_degree = c(2, 1))
#' print(resampled$resampled_data)
#' print(resampled$resampled_data_record)
#'
#' # Example 3: Resample data frame only
#' resampled <- get_resampled_data(data, resample_size = 10, resample_only = TRUE)
#' print(resampled$resampled_data)
#' print(resampled$resampled_data_record)
get_resampled_data <- function(data,
                               resample_size,
                               data_degree = NULL,
                               resample_only = FALSE) {
  if (is.null(data_degree)) data_degree <- rep(1, ncol(data))

  ## Resmapling Start
  # Resampling `data` (by columns) without any resampling method
  resampled_data <- as.data.frame(apply(data,
    2,
    sample,
    size = resample_size,
    replace = TRUE
  ))

  # Return resampled_data if resample_only == TRUE
  if (resample_only) {
    return(list(
      resampled_data = resampled_data,
      resampled_data_record = NULL
    ))
  }

  # Create data frame for recording the process of resampling each column of `data`
  resampled_data_record <- data.frame(
    Type = character(), # Record the type of each column of `data`
    Process = character() # Record the resampling process for each column of `data`
  )
  # TODO try apply/lapply
  for (n in 1:ncol(resampled_data)) {
    process <- "Coordinate"
    # Check continue or category
    if (is.continuous(data[, n])) {
      data_type <- "Continuous"
      # Check order and diversity
      if (data_degree[n] >= length(unique(data[, n])) - 1) {
        process <- c(process, "Smoothing")
        resampled_data[, n] <- do.call(
          get_resampling_methods,
          list(
            resampled_data[, n],
            data[, n],
            "smoothing"
          )
        )
      }
      # Check skewness
      if (abs(get_skewness(data[, n])) >= 1) {
        process <- c(process, "IQR")
        resampled_data[, n] <- do.call(
          get_resampling_methods,
          list(
            resampled_data[, n],
            data[, n],
            "IQR"
          )
        )
      }
    } else {
      data_type <- "Category"
      tbl <- table(data[, n])
      # Check data balance
      if (any(tbl < nrow(data) / (2 * length(tbl)))) {
        process <- c(process, "Flattening")
        resampled_data[, n] <- do.call(
          get_resampling_methods,
          list(
            resampled_data[, n],
            data[, n],
            "flattening"
          )
        )
      }
    }


    resampled_data_record <- rbind(resampled_data_record, data.frame(
      Type = data_type,
      Process = paste(process,
        collapse = "->"
      )
    ))
  }

  return(list(
    resampled_data = as.data.frame(resampled_data),
    resampled_data_record = as.data.frame(resampled_data_record)
  ))
}

#' Resampling Methods for Data
#'
#' This function applies various resampling methods to a vector `x` based on the chosen method.
#' It supports methods such as symmetrizing, smoothing, flattening, and IQR adjustment.
#'
#' @param x The vector to be resampled.
#' @param original_x The original vector before resampling.
#' @param method Character specifying the resampling method to apply. Options are "symmetrizing",
#'               "smoothing", "flattening", or "IQR".
#'
#' @return A vector of the same length as `x`, resampled according to the chosen method.
#' @export
#'
#' @examples
#' original_data <- c(10, 20, 10, 10, 20)
#' resampled_data <- get_resampling_methods(original_data, original_data, method = "smoothing")
#' print(resampled_data)
get_resampling_methods <- function(x,
                                   original_x,
                                   method = c(
                                     "symmetrizing",
                                     "smoothing",
                                     "flattening",
                                     "IQR"
                                   )) {
  method <- match.arg(method)
  return(
    switch(method,
      "symmetrizing" = mean(original_x) +
        sample(c(-1, 1), length(x), replace = TRUE) *
          (x - mean(original_x)),
      "smoothing" = x +
        stats::rnorm(length(x), 0, min(diff(sort(unique(original_x)))) / 10),
      "flattening" = (function(x) {
        x[sample(length(x), length(x) / 2)] <- sample(
          unique(original_x),
          length(x) / 2,
          replace = TRUE
        )
        return(x)
      })(x),
      "IQR" = (function(x) {
        x[sample(length(x), length(x) / 2)] <- rnorm(
          length(x) / 2,
          median(original_x),
          -get_iqr(original_x) / (2 * qnorm(1 / 4))
        )
        return(x)
      })(x)
    )
  )
}

#' Adjust `data_` Object from `data_initialization` Functions with Additional Matrix Calculations
#'
#' This function adjusts a data frame `data_` by adding new matrices based on the formula provided.
#' It computes adjusted matrices `adj_x`, `adj_obs_x`, and `adj_syn_x` and attaches them to the input `data_`.
#'
#' @param data_ A data frame containing the original data.
#' @param formula A formula specifying the adjustments to be made.
#'
#' @return The modified data frame `data_` with additional matrices `adj_x`, `adj_obs_x`, and `adj_syn_x`.
#' @export
get_adjusted_data <- function(data_, formula) {
  data_$adj_x <- model.frame(get_formula_rhs(as.formula(formula), TRUE), as.data.frame(data_$x))
  rownames(data_$adj_x) <- rownames(data_$x)
  data_$adj_obs_x <- as.data.frame(data_$adj_x[1:data_$obs_size, ])
  data_$adj_syn_x <- as.data.frame(data_$adj_x[(data_$obs_size + 1):data_$size, ])
  colnames(data_$adj_obs_x) <- colnames(data_$adj_syn_x) <- colnames(data_$adj_x)
  return(data_)
}

# ----- Functions for cat -----

#' Memoized Function to Retrieve Stan Models
#'
#' This function memoizes the creation of a Stan model based on the specified model type
#' and optional parameters for specific methods like GLMs. It returns a compiled Stan model object.
#'
#' @param model_type A character string specifying the type of Stan model (e.g., "glm", "lm").
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
get_stan_model <- memoise::memoise(function(model_type,
                                            glm_family = NULL,
                                            joint_tau = FALSE,
                                            glm_binomial_joint_theta = FALSE,
                                            glm_binomial_joint_alpha = FALSE,
                                            glm_binomial_gibbs = FALSE) {
  fam <- if (!is.null(glm_family)) get_glm_family(glm_family) else NULL

  return(rstan::stan_model(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      paste0(
        c(
          model_type,
          fam,
          if (any(
            joint_tau,
            glm_binomial_joint_theta,
            glm_binomial_joint_alpha
          )) {
            "joint"
          },
          if (model_type == "glm" && fam == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs) "theta",
          if (model_type == "glm" && fam == "binomial" && glm_binomial_joint_theta && !glm_binomial_gibbs & glm_binomial_joint_alpha) "alpha"
        ),
        collapse = "_"
      ),
      ".stan"
    )
  ))
})


# ----- Functions for cat_lmm -----

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
#' @export
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
update_sigma_square_cat <- function(sigma_square,
                                    sigma_u_square,
                                    d_seq,
                                    d_seq_star,
                                    q,
                                    q_star,
                                    tau) {
  return(sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) +
    tau / length(d_seq_star) * sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square)))
}


# ----- Functions for cat_cox -----

#' Compute the Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function calculates the Hessian matrix for a Cox proportional hazards model,
#' which is used to estimate the second-order partial derivatives of the log-likelihood function.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-event or censoring time for each observation.
#' @param event A numeric vector indicating whether the event of interest occurred (1) or the observation was censored (0).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector of entry times for each observation. Default is NULL, in which case all entry times are set to 0.
#'
#' @return A numeric matrix representing the Hessian matrix.
#' @export
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' event <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
#' coefs <- rnorm(ncol(X))
#'
#' # Compute Hessian matrix
#' hessian_matrix <- get_cox_hessian(X, time, event, coefs)
#' print(hessian_matrix)
get_cox_hessian <- function(X,
                            time,
                            event,
                            coefs,
                            entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  cersored_event <- which(as.numeric(event) == 1)
  hessian <- matrix(0, ncol(X), ncol(X))

  for (i in cersored_event) {
    risk_set_idx <- get_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      event_vector = event
    )
    risk_set <- matrix(X[risk_set_idx, ], ncol = ncol(X))

    exp_risk_lp <- exp(risk_set %*% coefs)
    risk_exp_risk_lp <- t(risk_set) %*% exp_risk_lp
    sum_exp_risk_lp <- sum(exp_risk_lp)

    risk_sqrt_exp_risk_lp <- sapply(
      1:length(risk_set_idx),
      function(i) (t(risk_set)[, i] * sqrt(exp_risk_lp[i]))
    )

    hessian <- hessian +
      crossprod(t(risk_sqrt_exp_risk_lp)) / sum_exp_risk_lp -
      tcrossprod(risk_exp_risk_lp) / sum_exp_risk_lp^2
  }

  return(-hessian)
}

#' Compute the Synthetic Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function calculates the synthetic Hessian matrix for a Cox proportional hazards model,
#' using the provided predictor variables, time-to-event data, coefficients, and a hazard constant.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-event or censoring time for each observation.
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param hazard_constant A numeric constant representing the hazard.
#'
#' @return A numeric matrix representing the synthetic Hessian matrix.
#' @export
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' coefs <- rnorm(ncol(X))
#' hazard_constant <- 0.1
#'
#' # Compute synthetic Hessian matrix
#' hessian_matrix <- get_syn_cox_hessian(X, time, coefs, hazard_constant)
#' print(hessian_matrix)
get_syn_cox_hessian <- function(X,
                                time,
                                coefs,
                                hazard_constant) {
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
#' time-to-event data, event indicators, coefficients, and optional entry points.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-event or censoring time for each observation.
#' @param event A numeric vector of event indicators (1 if event occurred, 0 if censored).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector of entry points for each observation. Default is NULL.
#'
#' @return A numeric vector representing the gradient for the Cox proportional hazards model.
#' @export
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' event <- c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1)
#' coefs <- rnorm(ncol(X))
#'
#' # Compute gradient
#' gradient <- get_cox_gradient(X, time, event, coefs)
#' print(gradient)
get_cox_gradient <- function(X,
                             time,
                             event,
                             coefs,
                             entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  gradient <- apply(
    cbind(time, X), 1,
    function(df) {
      df <- matrix(df, nrow = 1)
      ts <- df[, 1]
      xs <- df[, 2:ncol(df)]
      X_risk_set <- X[get_risk_set_idx(ts, entry_points, time, event), ] %>%
        matrix(ncol = ncol(X))

      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% coefs)
      t2 <- sum(exp(X_risk_set %*% coefs))

      return(xs - t1 / t2)
    }
  ) %*% event

  return(as.vector(gradient))
}

#' Compute the Gradient for Synthetic Data in a Cox Proportional Hazards Model
#'
#' This function calculates the gradient for synthetic data in a Cox proportional hazards model
#' using the provided predictor variables, time-to-event data, coefficients, and a hazard constant.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-event or censoring time for each observation.
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param hazard_constant A numeric value representing the hazard constant.
#'
#' @return A numeric vector representing the gradient for the synthetic data in the Cox proportional hazards model.
#' @export
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' coefs <- rnorm(ncol(X))
#' hazard_constant <- 0.01
#'
#' # Compute gradient
#' gradient <- get_syn_cox_gradient(X, time, coefs, hazard_constant)
#' print(gradient)
get_syn_cox_gradient <- function(X,
                                 time,
                                 coefs,
                                 hazard_constant) { # TODO Optimize code, e.g. replace for loop
  syn_lp <- as.numeric(get_linear_predictor(X, coefs))
  gradient <- rep(0, ncol(X))
  for (i in 1:nrow(X)) {
    gradient <- gradient + (1 - hazard_constant * time[i] * exp(syn_lp[i])) * X[i, ]
  }
  return(as.vector(gradient))
}

#' Compute the Partial Likelihood for a Cox Proportional Hazards Model
#'
#' This function calculates the partial likelihood for a Cox proportional hazards model
#' using the provided predictor variables, time-to-event data, event indicators, and coefficients.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param time A numeric vector representing the time-to-event or censoring time for each observation.
#' @param event A numeric vector indicating whether an event occurred (1) or the data is censored (0).
#' @param coefs A numeric vector of coefficients for the predictor variables.
#' @param entry_points An optional numeric vector representing the entry times for each observation. Default is NULL.
#'
#' @return A numeric value representing the partial likelihood for the Cox proportional hazards model.
#' @export
#'
#' @examples
#' # Example data
#' X <- matrix(rnorm(50), nrow = 10)
#' time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
#' event <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 1)
#' coefs <- rnorm(ncol(X))
#' entry_points <- NULL
#'
#' # Compute partial likelihood
#' pl <- get_cox_partial_likelihood(X, time, event, coefs, entry_points)
#' print(pl)
get_cox_partial_likelihood <- function(X,
                                       time,
                                       event,
                                       coefs,
                                       entry_points = NULL) {
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))

  cersored_event <- which(as.numeric(event) == 1)
  pl <- 0
  for (i in cersored_event) {
    risk_set_idx <- get_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      event_vector = event
    )
    risk_set <- matrix(X[risk_set_idx, ], ncol = ncol(X))
    exp_risk_lp <- exp(as.vector(get_linear_predictor(risk_set, coefs)))
    pl <- pl + sum(X[i, ] * coefs) - log(sum(exp_risk_lp))
  }
  return(pl)
}


get_cox_kappa <- function(X,
                          time,
                          event,
                          hazard_constant = NULL,
                          tol = 1e-5,
                          max_iter = 25,
                          entry_points = NULL,
                          initial_values = NULL) {
  if (is.null(hazard_constant)) hazard_constant <- sum(event) / mean(time) / nrow(X)
  if (is.null(entry_points)) entry_points <- rep(0, nrow(X))
  coefs <- if (is.null(initial_values)) rep(0, ncol(X)) else initial_values

  iter <- 1
  diff_new_old <- TRUE
  while ((iter < max_iter) & diff_new_old) {
    coefs_old <- coefs

    exp_lp <- exp(multiply_matrix(X, coefs_old))
    gd_old <- get_syn_cox_gradient(
      X = X,
      time = time,
      coefs = coefs_old,
      hazard_constant = hazard_constant
    )
    hess_old <- get_syn_cox_hessian(
      X = X,
      time = time,
      coefs = coefs_old,
      hazard_constant = hazard_constant
    )
    coefs <- coefs_old - get_qr_solve(hess_old, gd_old)

    iter <- iter + 1

    .true_stop(
      any(is.na(coefs)),
      "Please note NA appears in iteration."
    )

    .true_stop(
      norm(as.numeric(coefs), "2") > 1000000,
      "Please not the euclidean norm of coefficients is too large."
    )

    diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > tol * norm(as.numeric(coefs_old), "2")
  }

  lp_new <- as.numeric(get_linear_predictor(X, coefs))

  return(sum(lp_new - exp(lp_new) * hazard_constant * time) / nrow(X))
}

#' Compute Risk and Failure Sets for Survival Analysis
#'
#' This function computes the risk and failure sets for a given time vector, event vector, and a sequence of times (s).
#' It identifies the times at which events occur or subjects are censored.
#'
#' @param time_vector A numeric vector representing the time-to-event or censoring time for each observation.
#' @param event_vector A numeric vector indicating whether an event occurred (1) or the data is censored (0).
#' @param s A numeric vector of times at which to evaluate the risk and failure sets.
#'
#' @return A list containing two matrices: `risk_set` and `failure_set`, representing the risk and failure sets at each time point.
#' @export
#'
#' @examples
#' time_vector <- c(5, 8, 12, 18, 22)
#' event_vector <- c(1, 0, 1, 1, 0)
#' indicator_vector <- seq(1, 25, by = 1)
#' result <- get_risk_and_failure_sets(time_vector, event_vector, s)
#' print(result$risk_set)
#' print(result$failure_set)
get_risk_and_failure_sets <- function(time_vector,
                                      event_vector,
                                      indicator_vector) {
  # Function to find the time where the event occurred or reach censoring time
  find_occur_idx <- function(time) {
    min(which(s >= time))
  }

  risk_set <- failure_set <- matrix(0, length(time_vector), length(s))

  # Process subjects with event
  for (i in which(event_vector == 1)) {
    # Mark only exactly when the event occurrence
    failure_set[i, find_occur_idx(time_vector[i])] <- 1
    # Mark from risk up to the event time
    risk_set[i, 1:find_occur_idx(time_vector[i])] <- 1
  }

  # Process subjects without event
  for (i in which(event_vector == 0)) {
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
#' @param time_vector A numeric vector representing the time-to-event or censoring time for each subject.
#' @param event_vector A numeric vector indicating whether an event occurred (1) or the data is censored (0).
#'
#' @return A numeric vector of indices representing the subjects at risk at the specified time of interest.
#' @export
#'
#' @examples
#' time_of_interest <- 10
#' entry_vector <- c(0, 0, 0, 5, 7)
#' time_vector <- c(5, 8, 12, 18, 22)
#' event_vector <- c(1, 0, 1, 1, 0)
#' risk_set_idx <- get_risk_set_idx(time_of_interest, entry_vector, time_vector, event_vector)
#' print(risk_set_idx)
get_risk_set_idx <- function(time_of_interest,
                             entry_vector,
                             time_vector,
                             event_vector) {
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & event_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}

#' Solve a linear system using QR decomposition
#'
#' This function solves a linear system of equations represented by `coef_matrix %*% x = constant_vector`
#' using QR decomposition. It handles NA values in the solution by replacing them with a small constant.
#'
#' @param coef_matrix A matrix of coefficients from the linear system.
#' @param constant_vector The constant vector from the linear system.
#'
#' @return The solution vector `x` from the linear system.
#' @export
#'
#' @examples
#' coef_matrix <- matrix(c(1, 1, 1, 2, 3, 5, 4, 6, 8), nrow = 3, byrow = TRUE)
#' constant_vector <- c(10, 20, 30)
#' solution <- get_qr_solve(coef_matrix, constant_vector)
#' print(solution)
get_qr_solve <- function(coef_matrix,
                         constant_vector) {
  # Calculate solution using QR decomposition
  qr_solution <- solve(coef_matrix, constant_vector)

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
}

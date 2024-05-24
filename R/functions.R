#' Auto resampling of data
#'
#' This function performs automatic resampling on a dataset,
#' potentially applying different resampling methods based on
#' the type and characteristics of the data.
#'
#' @param data A matrix or data.frame, which is the original data for resampling.
#' @param resample_size A integer, which is the size(number of rows) of resampled data.
#' @param data_degree A numeric matrix or data.frame, which is the degree matrix of `data`.
#' @param resample_only A bool, which indicate if only do simple resmapling `data` or
#'                      resampling `data` with different resampling methods if need.
#' @param ...
#'
#' @return A list, including the resampled data in matrix format and
#'        the resampling process for future review
#' @export
#'
#' @examples
auto_resampling <- function(data,
                            resample_size,
                            data_degree = NULL,
                            resample_only = FALSE) {
  data <- as.matrix(data)
  resample_size <- as.integer(resample_size)
  data_degree <- if (is.null(data_degree)) rep(1, ncol(data)) else data_degree

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

  for (n in 1:ncol(resampled_data)) {
    process <- "Coordinate"
    # Check continue or category
    if (is.continuous(data[, n])) {
      data_type <- "Continuous"
      # Check order and diversity
      if (data_degree[n] >= length(unique(data[, n])) - 1) {
        process <- c(process, "Smoothing")
        resampled_data[, n] <- do.call(
          resampling_methods,
          list(
            resampled_data[, n],
            data[, n],
            "smoothing"
          )
        )
      }
      # Check skewness
      if (abs(e1071::skewness(data[, n])) >= 1) {
        process <- c(process, "IQR")
        resampled_data[, n] <- do.call(
          resampling_methods,
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
          resampling_methods,
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
    resampled_data = resampled_data,
    resampled_data_record = resampled_data_record
  ))
}



#' Resampling methods for one vector
#'
#'
#' This function applies various resampling methods
#' (symmetrizing, smoothing, flattening, and IQR adjustment)
#' to a vector `x` baesd on information from a vector `sample_x`.
#'
#' @param x A vector, which is the target to be resampled.
#' @param sample_x A vector, which is providing information (e.g. mean) for resmapling `x`.
#' @param method A character string, for specifying the resampling method to use.
#'
#' @return A vector, which is the resampled `x`.
#' @export
#'
#' @examples
#' # High skewed vector
#' x <- c(10000, rep(0, 9))
#' e1071::skewness(x)
#'
#' # Set sample x same as x in this example
#' sample_x <- x
#'
#' # Using one of resample method to reduce the skewness of x
#' resampled_x <- resampling_methods(x = x, sample_x = sample_x, method = "symmetrizing")
#' e1071::skewness(resampled_x)
#'
resampling_methods <- function(x,
                               sample_x,
                               method = c(
                                 "symmetrizing",
                                 "smoothing",
                                 "flattening",
                                 "IQR"
                               )) {
  method <- match.arg(method)
  return(
    switch(method,
      "symmetrizing" = mean(sample_x) +
        sample(c(-1, 1), length(x), replace = TRUE) *
          (x - mean(sample_x)),
      "smoothing" = x +
        stats::rnorm(length(x), 0, min(diff(sort(unique(sample_x)))) / 10),
      "flattening" = (function(x) {
        x[sample(length(x), length(x) / 2)] <- sample(
          unique(sample_x),
          length(x) / 2,
          replace = TRUE
        )
        return(x)
      })(x),
      "IQR" = (function(x) {
        x[sample(length(x), length(x) / 2)] <- rnorm(
          length(x) / 2,
          median(sample_x),
          -get_iqr_diff(sample_x) / (2 * qnorm(1 / 4))
        )
        return(x)
      })(x)
    )
  )
}

get_log_density <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          Y,
                          model = NULL,
                          weights = 1,
                          ...,
                          coefficients = NULL) {
      return(log_density_binomial(
        Y = as.matrix(Y),
        mean = multiply_two_matrices(
          X,
          if (is.null(model)
          ) {
            coefficients
          } else {
            model$coefficients
          }
        ),
        weights = weights
      ))
    },
    "gaussian" = function(model = NULL,
                          weights = 1,
                          ...,
                          residuals = NULL) {
      return(log_density_gaussian(
        residuals = if (is.null(model)
        ) {
          residuals
        } else {
          model$residuals
        },
        weights = weights
      ))
    }
  )
}


get_diag_approx_cov <- function(X, model) {
  return(diag(chol2inv(chol(second_derivative(
    X = as.matrix(X),
    fitted_values = model$fitted.values,
    prior_weights = model$prior.weights
  )))))
}


get_log_density_grad <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          Y,
                          model = NULL,
                          weights = 1,
                          ...,
                          coefficients = NULL) {
      return(
        log_density_grad_binomial(
          X = as.matrix(X),
          Y = as.matrix(Y),
          mean = get_mean("binomial")(
            X = X,
            coefficients = if (is.null(model)
            ) {
              coefficients
            } else {
              model$coefficients
            }),
          weights = weights
        )
      )
    }
  )
}

get_mean <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          model = NULL,
                          ...,
                          coefficients = NULL) {
      return(
        mean_of_binomial(
          X = as.matrix(X),
          coefficients = if (is.null(coefficients)
          ) {
            coefficients
          } else {
            coef(model)
          }
        )
      )
    },
    "gaussian" = function(X,
                          model = NULL,
                          ...,
                          coefficients = NULL) {
      return(
        mean_of_gaussian(
          X = as.matrix(X),
          coefficients = as.matrix(if (is.null(coefficients)
          ) {
            coef(model)
          } else {
            coefficients
          })
        )
      )
    }
  )
}


get_sim_y <- function(fam) {
  switch(fam,
    "binomial" = function(X,
                          model = NULL,
                          ...,
                          coefficients = NULL,
                          mean = NULL) {
      get_data(
        "binomial"
      )(
        n = nrow(X),
        mean = if (is.null(mean)) {
          get_mean(
            "binomial"
          )(
            X = X,
            model = model,
            coefficients = coefficients
          )
        } else {
          mean
        })
    },
    "gaussian" = function(X,
                          model = NULL,
                          coefficients = NULL,
                          ...,
                          mean = NULL,
                          sd = NULL) {
      get_data(
        "gaussian"
      )(
        n = nrow(X),
        mean = if (is.null(mean)) {
          get_mean(
            "gaussian"
          )(
            X = X,
            coefficients = coefficients,
            model = model
          )
        } else {
          mean
        },
        sd = if (is.null(sd)) sqrt(sum(residuals(model)^2) / df.residual(model)) else sd
      )
    }
  )
}

get_lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, coefficients, fam, ...) 2 * get_mean(fam)(X, coefficients),
    "classification_error" = function(X, coefficients, fam, ...) 2 * get_mean(fam)(X, coefficients),
    "logistic_deviance" = function(X, coefficients, fam, ...) {
      .not_same_stop(
        fam,
        "binomial",
        "only binomial family can use logistic_deviance"
      )
      return(multiply_two_matrices(X, coefficients))
    }
  )
}

get_iqr_diff <- function(x, invl = 0.5) {
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


get_data <- function(fam) {
  switch(fam,
    "binomial" = function(n = 10,
                          mean = 0,
                          ...) {
      stats::rbinom(
        n,
        1,
        mean
      )
    },
    "gaussian" = function(n = 10,
                          mean = 0,
                          sd = 1,
                          ...) {
      stats::rnorm(
        n,
        mean,
        sd
      )
    },
    "uniform" = function(n = 10,
                         mean = 0,
                         ...) {
      rep(
        mean,
        n
      )
    }
  )
}

get_discrepancy <- function(func) {
  switch(func,
    "logarithmic_error" = function(Y,
                                   est_Y = NULL,
                                   ...,
                                   fam = NULL,
                                   est_model = NULL,
                                   est_coefficients = NULL) {
      # cap the value between 0.0001 and 0.9999 for better calculation process
      browser()
      return(
        discrepancy_of_logarithmic_error(
          Y = pmax(0.0001, pmin(0.9999, Y)),
          est_Y = pmax(0.0001, pmin(
            0.9999,
            if (is.null(est_Y)) {
              get_mean(
                fam
              )(
                X = X,
                coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
              )
            } else {
              est_Y
            }
          ))
        )
      )
    },
    "square_error" = function(Y,
                              est_Y = NULL,
                              ...,
                              X = NULL,
                              fam = "gaussian",
                              est_model = NULL,
                              est_coefficients = NULL) {
      return(
        discrepancy_of_square_error(
          Y = Y,
          est_Y = if (is.null(est_Y)) {
            get_mean(
              fam
            )(
              X = X,
              coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
            )
          } else {
            est_Y
          }
        )
      )
    },
    "classification_error" = function(Y,
                                      est_Y = NULL,
                                      ...,
                                      X = NULL,
                                      fam = NULL,
                                      est_model = NULL,
                                      est_coefficients = NULL) {
      return(
        discrepancy_of_classification_error(
          Y = Y,
          est_Y = if (is.null(est_Y)) {
            get_mean(
              fam
            )(
              X = X,
              coefficients = if (is.null(est_model)) est_coefficients else est_model$coefficients
            )
          } else {
            est_Y
          }
        )
      )
    },
    "logistic_deviance" = function(X,
                                   Y,
                                   fam,
                                   model = NULL,
                                   coefficients = NULL,
                                   ...) {
      -get_log_density(fam)(
        X = X,
        Y = Y,
        model = model,
        coefficients = coefficients)
    }
  )
}

update_sigma_square <- function(sigma_square,
                                sigma_u_square,
                                d_seq,
                                q) { # minimization program
  sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
}


update_sigma_u_square <- function(sigma_u_square,
                                  sigma_square,
                                  d_seq,
                                  q) { # minimization program
  sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
}


update_sigma_square_cat <- function(sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau) { # minimization program
  sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) + tau / length(d_seq_star) *
    sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))
}

update_sigma_u_square_cat <- function(sigma_u_square, sigma_square, d_seq, d_seq_star, q, q_star, tau) { # minimization program
  sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) + tau / length(d_seq_star) *
    sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))
}


stan_model_selection <- memoise::memoise(function(fam,
                                                  joint,
                                                  gibbs,
                                                  use_theta,
                                                  joint_alpha = FALSE) {
  return(rstan::stan_model(
    file = if (is.null(model_file)) {
      paste0(
        "~/catalytic/R/rstan_script/",
        fam,
        if (joint & isFALSE(gibbs)) "_joint",
        "_bayesian",
        if (fam == "binomial" & joint & use_theta & isFALSE(gibbs)) "_theta",
        if (fam == "binomial" & joint & isFALSE(use_theta) & isFALSE(gibbs)) "_tau",
        if (fam == "binomial" & joint & use_theta & isFALSE(gibbs) & joint_alpha) "_joint_alpha",
        ".stan"
      )
    } else {
      model_file
    }
  ))
})

standardize_data <- function(data, if_missing = NULL) {
  data <- as.data.frame(data)

  for (i in 1:ncol(data)) {
    if (!(NA %in% suppressWarnings(as.numeric(data[, i])))) {
      data[, i] <- as.numeric(data[, i])
    } else if (!(NA %in% suppressWarnings(as.logical(data[, i])))) {
      data[, i] <- as.logical(data[, i])
    } else {
      data[, i] <- as.factor(data[, i])
    }
  }
  return(data)
}


validate_glm_formula <- function(f) {
  if (!any(grepl("~", f))) {
    stop("please enter ~ for formula",
      call. = FALSE
    )
  }
}

validate_surv_formula <- function(f) {
  validate_glm_formula(f)

  if (!any(grepl("Surv", f))) {
    stop("please use `Surv` for formula",
      call. = FALSE
    )
  }
}



validate_family <- function(f) {
  if (is.character(f)) {
    f <- get(f, mode = "function", envir = parent.frame(2))
  }
  if (is.function(f)) {
    f <- f()
  }
  if (!is(f, "family")) {
    stop("'family' must be a family.", call. = FALSE)
  }

  return(f)
}

get_family <- function(family) {
  fam <- try(family$family, silent = TRUE)
  if (inherits(fam, "try-error")) {
    fam <- family()$family
  }
  return(fam)
}

family_plus_link <- function(fam) {
  if (is.character(fam)) {
    stopifnot(identical(fam, x$method))
    fam <- paste0("ordered [", fam, "]")
  } else {
    fam <- paste0(fam$family, " [", fam$link, "]")
  }
  return(fam)
}

formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) {
    "--MARK--"
  } else {
    " "
  }
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent) {
    return(char)
  }
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

get_rhs <- function(formula, with_ = FALSE) {
  # Split the formula string by "~" and extract the first element
  if (with_) {
    return(as.formula(paste("~", deparse(rlang::f_rhs(formula)), collapse = " ")))
  }
  return(rlang::f_rhs(formula))
}

get_lhs <- function(formula) {
  return(rlang::f_lhs(formula))
}

H_matrix <- function(X, X_star, tau, M, n) {
  X <- as.matrix(X)
  X_star <- as.matrix(X_star)
  return(X %*% (
    chol2inv(chol((t(X) %*% X + (tau / M) * (t(X_star) %*% X_star))))) %*% (
    t(X) + tau / (M * n) * t(X_star) %*% rep(1, M) %*% t(rep(1, n))
  ))
}


get_adj_data <- function(data_, formula) {
  data_$adj_x <- model.frame(get_rhs(formula, TRUE), as.data.frame(data_$x))
  data_$adj_obs_x <- data_$adj_x[1:data_$obs_size, ]
  data_$adj_syn_x <- data_$adj_x[(data_$obs_size + 1):data_$size, ]
  return(data_)
}



check_incul_wo_overlap <- function(whole_set, ...) {
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


#' Solves a system of linear equations using QR decomposition.
#'
#' This function solves the system of linear equations \(Mx = b\) using the QR decomposition method.
#' If any NA values are encountered in the solution, they are replaced with 0.0001.
#'
#' @param coef_matrix A matrix of coefficients.
#' @param constant_vector A vector representing the constants on the right side of the equations.
#'
#' @return A vector containing the solutions to the system of equations.
#'
#'
#' @seealso \code{\link[base]{qr}} for the QR decomposition function used internally.
get_qr_solve <- function(coef_matrix, constant_vector) {
  # Calculate solution using QR decomposition
  qr_solution <- qr.coef(qr(coef_matrix), constant_vector)

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
}


# Hessian function --------------------------------------------------------
#' Compute the Hessian Matrix for Log Partial Likelihood of Cox Proportional Hazards Model
#'
#' This function calculates the negative Hessian matrix associated with the log partial likelihood of the Cox proportional hazards model.
#'
#' @param beta A vector representing the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param entry A vector of entry times for each observation.
#' @param Ts A vector of event times for each observation.
#' @param event A binary vector indicating event occurrence (1 for event, 0 for censored).
#'
#' @return The negative Hessian matrix for the log partial likelihood of the Cox proportional hazards model.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
get_cox_hessian <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  hessian <- matrix(0, p, p)
  uncensored_set <- which(as.numeric(event) == 1)

  for (i in uncensored_set) {
    ts <- Ts[i]
    xs <- Xs[i, ]
    riskindex <- get_risk_set(ts, entry, Ts, event)

    if (length(riskindex) > 1) {
      X_risk_set <- Xs[riskindex, ]
    } else {
      X_risk_set <- matrix(Xs[riskindex, ], ncol = ncol(Xs))
    }

    theta <- exp(X_risk_set %*% beta)
    t1 <- t(X_risk_set) %*% theta
    t2 <- sum(theta)

    xxxt <- t(X_risk_set)
    eee <- sapply(1:length(riskindex), function(i) (xxxt[, i] * sqrt(theta[i])))

    hessian <- hessian + crossprod(t(eee)) / t2 - t1 %*% t(t1) / t2^2
  }

  return(-hessian)
}


#' Calculate the Hessian Matrix for the Log Catalytic with No Censoring for Synthetic Data
#'
#' This function computes the Hessian matrix for the log catalytic prior, specifically for synthetic data that has no censoring.
#'
#' @param Xsbeta A vector representing the product of the covariates and the coefficients of the model.
#' @param Xs A matrix of covariates (predictor variables) where rows represent observations and columns represent features.
#' @param Ts A vector of event times for each observation.
#' @param h_daga A function representing the hazard rate for synthetic data.
#'
#' @return A Hessian matrix for the log catalytic prior.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
get_syn_cox_hessian <- function(Xsbeta, Xs, Ts, h_daga) {
  M <- nrow(Xs)
  p <- ncol(Xs)
  hess <- matrix(0, p, p)

  for (ii in 1:M) {
    hess <- hess - integrate(Vectorize(h_daga), lower = 0, upper = Ts[ii])$value * exp(Xsbeta[ii]) * tcrossprod(Xs[ii, ])
  }

  return(hess)
}


get_risk_set <- function(time_of_interest,
                         entry_vector,
                         time_vector,
                         event_vector) {
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & event_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}



# gradient function -------------------------------------------------------
#' Compute the Gradient for Log Partial Likelihood of Cox Proportional Hazards Model
#'
#' This function calculates the gradient vector associated with the log partial likelihood of the Cox proportional hazards model.
#'
#' @param beta A vector representing the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param entry A vector of entry times for each observation.
#' @param Ts A vector of event times for each observation.
#' @param event A binary vector indicating event occurrence (1 for event, 0 for censored).
#'
#' @return A vector representing the gradient of the log partial likelihood with respect to beta.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
get_cox_gradient <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)

  gradient <- apply(cbind(Ts, Xs), 1, function(df) {
    df <- matrix(df, nrow = 1)
    ts <- df[, 1]
    xs <- df[, 2:ncol(df)]
    X_risk_set <- Xs[get_risk_set(ts, entry, Ts, event), ] %>%
      matrix(ncol = ncol(Xs))

    t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
    t2 <- sum(exp(X_risk_set %*% beta))

    return(xs - t1 / t2)
  }) %*% event

  return(gradient)
}



#' Compute the Gradient for Log Catalytic Prior for Synthetic Data
#'
#' This function calculates the gradient of the log catalytic prior, specifically for synthetic data without censoring.
#'
#' @param Xsbeta A vector representing the product of the covariates and the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param Ts A vector of event times for each observation.
#' @param status.s A vector indicating the status of the synthetic data (typically 1 for event and 0 for censored).
#' @param h_daga A function representing the hazard rate for synthetic data.
#'
#' @return A vector representing the gradient of the log catalytic prior with respect to beta.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
get_syn_cox_gradient <- function(Xsbeta, Xs, Ts, status.s, h_daga) {
  M <- nrow(Xs)
  p <- ncol(Xs)
  gd <- rep(0, p)

  for (ii in 1:M) {
    gd <- gd + (status.s[ii] - integrate(Vectorize(h_daga), lower = 0, upper = Ts[ii])$value * exp(Xsbeta[ii])) * Xs[ii, ]
  }

  return(gd)
}


setting.interval <- function(y, delta, s, J) {
  n <- length(y)

  smax <- max(s)

  case0 <- which(delta == 0)
  case1 <- which(delta == 1)

  case0yleq <- which(delta == 0 & y <= smax)
  case0ygeq <- which(delta == 0 & y > smax)
  case1yleq <- which(delta == 1 & y <= smax)
  case1ygeq <- which(delta == 1 & y > smax)


  ind.d <- ind.r <- matrix(0, n, J)

  for (i in case1yleq) {
    d.mat.ind <- min(which(s - y[i] >= 0))
    ind.d[i, d.mat.ind] <- 1
    ind.r[i, 1:d.mat.ind] <- 1
  }

  for (i in case0yleq) {
    cen.j <- min(which(s - y[i] >= 0))
    ind.r[i, 1:cen.j] <- 1
  }

  if (length(union(case1ygeq, case0ygeq)) > 0) {
    ind.r[union(case1ygeq, case0ygeq), ] <- 1
  }

  ind.r_d <- ind.r - ind.d


  d <- colSums(ind.d)

  list(
    ind.r = ind.r,
    ind.d = ind.d,
    d = d,
    ind.r_d = ind.r_d
  )
}



real_coefficients_vector <- function(dimention,
                                     zeta,
                                     r,
                                     pseudo_size = 1000) {
  # generate real beta vector shown in paper

  real_coefficients0 <- c(rep(0, dimention))

  # ζ zeta proportion of the coordinates of β(0) is randomly selected and set to 1
  real_coefficients0[sample(dimention, (dimention) * zeta, replace = FALSE)] <- 1
  real_coefficients0 <- c(0, rnorm(dimention) * real_coefficients0)
  pseudo_covariate <- observed_covariate_vector(dimention, pseudo_size, 0, diag(rep(1, dimention)), 0.5)
  coefficients_scale <- 1

  repeat {
    real_coefficients <- real_coefficients0 * coefficients_scale
    real_coefficients[1] <- -mean(pseudo_covariate %*% real_coefficients[-1])
    pseudo_response <- get_mean("binomial")(pseudo_covariate, coefficients = real_coefficients)
    oracle_error <- mean(1 - abs(1 - 2 * pseudo_response)) / 2
    if (abs(oracle_error - r) < 0.005) {
      break
    }
    coefficients_scale <- ifelse(oracle_error > r, coefficients_scale * 1.1, coefficients_scale * 0.9)
  }
  return(real_coefficients)
}

observed_covariate_vector <- function(dimention, # number of columns
                                      sample_size, # number of rows
                                      Z_mean,
                                      Z_variance,
                                      Z_corr) {
  # generate observed covariate vector from Gaussian random vector Z

  X_i <- Z_i <- MASS::mvrnorm(
    n = sample_size,
    mu = rep(Z_mean, dimention),
    Sigma = correlation_matrix(Z_corr, dimention) %*% Z_variance
  )


  # set Xi,j = 2 · 1Zi,j >0 − 1, 2j < p Zi,j , 2j ≥ p
  X_i[, 1:floor(dimention / 2)] <- 2 * (Z_i[, 1:floor(dimention / 2)] > 0) - 1

  return(X_i)
}

err <- function(x,
                coefficients0,
                coefficients1) {
  return(mean(get_discrepancy("logarithmic_error")(Y = get_mean("binomial")(x, coefficients = coefficients0), est_Y = get_mean("binomial")(x, coefficients = coefficients1))))
}

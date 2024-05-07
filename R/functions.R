source("~/catalytic/R/utils.R")

auto_resampling <- function(obs_data,
                            resample_size,
                            data_degree = NULL, ...) {
  recd_resampled_data <- data.frame(
    Type = character(),
    Process = character(),
    stringsAsFactors = FALSE
  )

  data_degree <- if (is.null(data_degree)) rep(1, ncol(obs_data)) else data_degree

  resampled_data <- as.data.frame(apply(obs_data,
    2,
    sample,
    size = resample_size,
    replace = TRUE
  ))

  for (n in 1:ncol(resampled_data)) {
    process <- "Coordinate"
    # Check continue or category
    if (is.continuous(obs_data[, n])) {
      data_type <- "Continuous"
      # Check order and diversity
      if (data_degree[n] >= length(unique(obs_data[, n])) - 1) {
        process <- c(process, "Smoothing")
        resampled_data[, n] <- do.call(
          resampling_methods,
          list(
            resampled_data[, n],
            obs_data[, n],
            "smoothing"
          )
        )
      }
      # Check skewness
      if (abs(e1071::skewness(obs_data[, n])) >= 1) {
        process <- c(process, "IQR")
        resampled_data[, n] <- do.call(
          resampling_methods,
          list(
            resampled_data[, n],
            obs_data[, n],
            "IQR"
          )
        )
      }
    } else {
      data_type <- "Category"
      tbl <- table(obs_data[, n])
      # Check data balance
      if (any(tbl < nrow(data) / (2 * length(tbl)))) {
        process <- c(process, "Flattening")
        resampled_data[, n] <- do.call(
          resampling_methods,
          list(
            resampled_data[, n],
            obs_data[, n],
            "flattening"
          )
        )
      }
    }


    recd_resampled_data <- rbind(recd_resampled_data, data.frame(
      Type = data_type,
      Process = paste(process,
        collapse = "->"
      )
    ))
  }

  return(list(
    resampled_data = resampled_data,
    recd_resampled_data = recd_resampled_data
  ))
}



resampling_methods <- function(x, sample_x, method) {
  switch(method,
    "symmetrizing" = mean(sample_x) + sample(c(-1, 1),
      length(x),
      replace = TRUE
    ) * (x - mean(sample_x)),
    "smoothing" = x + stats::rnorm(length(x), 0, min(diff(sort(unique(sample_x)))) / 10),
    "flattening" = (function(x) {
      x[sample(length(x), length(x) / 2)] <- sample(unique(sample_x),
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
          coefficients = if (is.null(coefficients)
          ) {
            coef(model)
          } else {
            coefficients
          }
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
                              fam = NULL,
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
  data_$adj_x <- model.frame(get_rhs(formula, TRUE), data_$x)
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

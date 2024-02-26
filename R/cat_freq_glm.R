#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param risk_estimate_method
#' @param discrepancy_function
#' @param taus
#' @param subset
#' @param na.action
#' @param offset
#' @param ...
#' @param tau_0
#' @param bootstrap_iteration_times
#' @param cv_fold_num
#'
#' @return
#' @export
#'
#' @examples
cat_freq_glm <-
  function(formula,
           family = gaussian(),
           data_,
           risk_estimate_method = c(
             "parametric_bootstrap",
             "cross_validation",
             "mallows_estimate",
             "steinian_estimate"
           ),
           discrepancy_function = c(
             "square_error",
             "classification_error",
             "logistic_deviance"
           ),
           taus = NULL,
           subset = NULL,
           na.action = NULL,
           offset = NULL,
           ...,
           tau_0 = 1,
           bootstrap_iteration_times = 100,
           cv_fold_num = 5) {
    .not_in_stop(
      c(
        "obs_size", "obs_data", "obs_x", "obs_y", "syn_size", "syn_data", "syn_x",
        "syn_y", "syn_model", "data", "size", "x", "y"
      ),
      names(data_),
      "Not all necessary variables are present in data_."
    )

    risk_estimate_method <- match.arg(risk_estimate_method)
    discrepancy_function <- match.arg(discrepancy_function)
    validate_glm_formula(formula)
    family <- validate_family(family)

    data_$adj_x <- model.frame(rhs(formula, TRUE), data_$x)
    data_$adj_obs_x <- data_$adj_x[1:data_$obs_size, ]
    data_$adj_syn_x <- data_$adj_x[(data_$obs_size + 1):data_$size, ]

    if (is.null(taus)) taus <- seq(max(0, ncol(data_$data) / 4 - 2), ncol(data_$data) / 4 + 2, 0.5)

    validate_tau(taus)

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "family", "subset", "na.action", "offset"),
      table = names(mf), nomatch = 0L
    )
    mf <- mf[c(1L, m)]


    .large_number_warning(
      data_$obs_size,
      if (risk_estimate_method == "parametric_bootstrap") 100 else if (
        risk_estimate_method == "steinian_estimate") 50 else Inf,
      "Number of observation data is large, suggest use cross-validation."
    )

    .large_number_ask(
      data_$obs_size,
      if (risk_estimate_method == "parametric_bootstrap") 150 else if (
        risk_estimate_method == "steinian_estimate") 100 else Inf,
      "Number of observation data is too large, do you want to continue? y/n",
      "n"
    )


    min_risk_estimate <- Inf
    risk_estimate_collect <- c()
    optimal_tau <- NULL

    for (tau in taus) {
      risk_estimate_collect <- append(
        risk_estimate_collect,
        do.call(
          risk_estimate_method,
          c(
            tau = tau,
            as.list(call[-1])
          )
        )
      )
    }

    optimal_tau <- taus[which.min(risk_estimate_collect)]

    mf$data <- data_$data
    mf$weight <- c(rep(1, data_$obs_size), rep(optimal_tau / data_$syn_size, data_$syn_size))

    freq_model <- do.call(stats::glm, as.list(mf[-1]))

    call["data_"] <- NULL
    call["taus"] <- NULL
    call["family"] <- NULL
    cat_freq <- c(list(
      data_ = data_,
      taus = taus,
      family = family,
      optimal_tau = optimal_tau,
      optimal_freq_model = freq_model,
      optimal_beta = freq_model$coefficients
    ), as.list(call[-1]))

    class(cat_freq) <- c(cat_freq$class, "cat_freq_glm")
    return(cat_freq)
  }


parametric_bootstrap <- function(formula,
                                 family = gaussian(),
                                 data_,
                                 tau,
                                 subset = NULL,
                                 na.action = NULL,
                                 offset = NULL,
                                 ...,
                                 discrepancy_function = "square_error",
                                 tau_0 = 1,
                                 bootstrap_iteration_times = 100) {
  fam <- family$family
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "family", "subset", "na.action", "offset"),
    table = names(mf), nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))
  mf$data <- data_$data
  beta_tau_model <- do.call(stats::glm, as.list(mf[-1]))


  D <- mean(discrepancy_functions(discrepancy_function)(
    data_$adj_obs_x,
    data_$obs_y,
    beta_tau_model$coefficients,
    fam)[[1]])

  mf$weight <- c(rep(1, data_$obs_size), rep(tau_0 / data_$syn_size, data_$syn_size))
  hat_tau_0_beta_model <- do.call(stats::glm, as.list(mf[-1]))

  mean_Y <- fam_mean(fam)(
    data_$adj_obs_x,
    hat_tau_0_beta_model$coefficients)

  boot_Y_collect <-
    lambda_boot_hat_mu_collect <- matrix(NA,
      nrow = data_$obs_size,
      ncol = bootstrap_iteration_times
    )

  for (bootstrap_iter in 1:bootstrap_iteration_times) {
    boot_Y <- generate_data_from_dist(
      fam,
      data_$obs_size,
      mean_Y
    )
    boot_Y_df <- as.data.frame(boot_Y)
    names(boot_Y_df) <- names(data_$syn_y)
    mf$data <- cbind(data_$x, rbind(boot_Y_df, data_$syn_y))
    mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))

    boot_hat_tau_beta_model <- do.call(stats::glm, as.list(mf[-1]))
    boot_Y_collect[, bootstrap_iter] <- boot_Y

    lambda_boot_hat_mu_collect[, bootstrap_iter] <- lambda_function(
      discrepancy_function
    )(data_$adj_obs_x,
      boot_hat_tau_beta_model$coefficients,
      fam)
  }

  return(D + sum(sweep(boot_Y_collect, 1, rowMeans(boot_Y_collect)) *
    lambda_boot_hat_mu_collect) / bootstrap_iteration_times)
}

steinian_estimate <- function(formula,
                              family = gaussian(),
                              data_,
                              tau,
                              subset = NULL,
                              na.action = NULL,
                              offset = NULL,
                              ...,
                              discrepancy_function = "square_error",
                              tau_0 = 1) {
  fam <- family$family
  .not_same_stop(
    fam,
    "binomial",
    "Only binomial can use steinian_estimate."
  )


  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "family", "subset", "na.action", "offset"),
    table = names(mf), nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))
  mf$data <- data_$data

  beta_tau_model <- do.call(stats::glm, as.list(mf[-1]))

  D <- mean(discrepancy_functions(discrepancy_function)(data_$obs_x,
    data_$adj_obs_y,
    beta_tau_model$coefficients,
    fam)[[1]])

  mf$weight <- c(rep(1, data_$obs_size), rep(tau_0 / data_$syn_size, data_$syn_size))
  hat_beta_0_model <- do.call(stats::glm, as.list(mf[-1]))

  hat_mu_0 <- fam_mean(fam)(data_$adj_obs_x, hat_beta_0_model$coefficients)

  with_ <- with(hat_beta_0_model, fitted.values * (1 - fitted.values))

  steinian_sum <- 0
  for (obs in 1:data_$obs_size) {
    triangle_response <- data_$y
    triangle_response[obs] <- 1 - data_$y[obs]

    mf$data <- cbind(data_$adj_x, triangle_response)
    mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))
    trigangle_hat_tau_beta_model <- do.call(stats::glm, as.list(mf[-1]))

    steinian_sum <- steinian_sum +
      with_[obs] * (2 * data_$y[obs] - 1) * (tcrossprod(
        as.matrix(data_$adj_x[obs, ]),
        as.matrix(beta_tau_model$coefficient - trigangle_hat_tau_beta_model$coefficient)[-1]
      ) + (
        beta_tau_model$coefficient - trigangle_hat_tau_beta_model$coefficient)[1])
  }

  return(D + steinian_sum)
}

mallows_estimate <- function(formula,
                             family = gaussian(),
                             data_,
                             tau,
                             subset = NULL,
                             na.action = NULL,
                             offset = NULL,
                             ...,
                             tau_0 = 1) {
  fam <- family$family
  .not_same_stop(
    fam,
    "gaussian",
    "Only gaussian can use mallows_estimate"
  )

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "family", "subset", "na.action", "offset"),
    table = names(mf), nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))
  mf$data <- data_$data

  beta_tau_model <- do.call(stats::glm, as.list(mf[-1]))

  D <- mean(discrepancy_functions("square_error")(data_$obs_x,
    data_$adj_obs_y,
    beta_tau_model$coefficient,
    fam)[[1]])

  mf$weight <- c(rep(1, data_$obs_size), rep(tau_0 / data_$syn_size, data_$syn_size))
  beta_0_model <- do.call(stats::glm, as.list(mf[-1]))

  var_Y <- sum((y_value(fam)(data_$adj_obs_x, beta_0_model$coefficient) - data_$obs_y)^2) / data_$obs_size

  tau_H <- H_matrix(
    X = data_$adj_obs_x,
    X_star = data_$adj_syn_x,
    tau = tau,
    M = data_$syn_size,
    n = data_$obs_size
  )

  return(D + 2 * var_Y * sum(diag(tau_H)) / data_$obs_szie)
}


H_matrix <- function(X, X_star, tau, M, n) {
  X <- as.matrix(X)
  X_star <- as.matrix(X_star)
  return(X %*% (
    (t(X) %*% X + (tau / M) * (t(X_star) %*% X_star))^(-1)) %*% (
    t(X) + tau / (M * n) * t(X_star) %*% rep(1, M) %*% t(rep(1, n))
  ))
}


cross_validation <- function(formula,
                             family = gaussian(),
                             data_,
                             tau,
                             subset = NULL,
                             na.action = NULL,
                             offset = NULL,
                             ...,
                             discrepancy_function = "square_error",
                             cv_fold_num = 5) {
  error_sum <- 0
  random_sample <- sample(1:data_$obs_size)

  for (fold in 1:cv_fold_num) {
    test_x <- data_$adj_obs_x[random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]
    test_y <- data_$obs_y[random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)]]

    train_x <- data_$adj_obs_x[-random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]
    train_y <- data_$obs_y[-random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)]]

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "family", "subset", "na.action", "offset"),
      table = names(mf), nomatch = 0L
    )
    mf <- mf[c(1L, m)]
    mf$weight <- c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size))
    mf$data <- rbind(cbind(train_x, train_y), data_$syn_data)
    estimate_beta_model <- do.call(stats::glm, as.list(mf[-1]))

    error_sum <- error_sum + mean(discrepancy_functions(discrepancy_function)(
      test_x,
      test_y,
      estimate_beta_model$coefficient,
      fam)[[1]])
  }

  return(error_sum / cv_fold_num)
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_freq_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nRisk estimate method: ", x$risk_estimate_method,
    "\n\nOptimal tau: ", x$optimal_tau,
    "\n\nOptimal estimated beta: ", x$optimal_beta
  )
  invisible(x)
}

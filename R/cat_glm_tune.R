#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param risk_estimate_method
#' @param discrepancy_function
#' @param taus
#' @param ...
#' @param tau_0
#' @param bootstrap_iteration_times
#' @param cv_fold_num
#'
#' @return
#' @export
#'
#' @examples
cat_glm_tune <-
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
           ...,
           tau_0 = NULL,
           bootstrap_iteration_times = 100,
           cv_fold_num = 5) {
    .not_true_stop(
      is.data_(data_),
      paste(
        "Please enter `data_` that generated from",
        "catalytic::data_initialization."
      )
    )

    validate_glm_formula(formula)
    family <- validate_family(family)
    fam <- get_family(family)

    if (missing(risk_estimate_method)) {
      risk_estimate_method <- if (
        data_$obs_size > 200) {
        "cross_validation"
      } else {
        "parametric_bootstrap"
      }
    }

    if (missing(discrepancy_function)) {
      discrepancy_function <- if (
        fam == "binomial") {
        "logistic_deviance"
      } else {
        "square_error"
      }
    }

    risk_estimate_method <- match.arg(risk_estimate_method)
    discrepancy_function <- match.arg(discrepancy_function)

    data_ <- get_adj_data(data_ = data_, formula = formula)

    if (is.null(taus)) taus <- seq(max(0, ncol(data_$data) / 4 - 2), ncol(data_$data) / 4 + 2, 0.5)
    if (is.null(tau_0)) tau_0 <- if (fam == "gaussian") 1 else ncol(data_$adj_x) / 4

    .small_number_stop(taus, 0, "Please enter non-negative `taus` value.")
    .small_number_stop(tau_0, 0, "Please enter non-negative `tau_0` value.")

    .large_number_warning(
      data_$obs_size,
      if (risk_estimate_method == "parametric_bootstrap") {
        150
      } else if (
        risk_estimate_method == "steinian_estimate") {
        100
      } else {
        Inf
      },
      paste(
        "Due to the considerable size of the dataset,",
        "Please consider utilizing `risk_estimate_method = cross_validation`",
        "for quicker estimate process."
      )
    )

    .large_number_ask(
      data_$obs_size,
      if (risk_estimate_method == "parametric_bootstrap") {
        200
      } else if (
        risk_estimate_method == "steinian_estimate") {
        150
      } else {
        Inf
      },
      paste0(
        "The large size of the dataset ",
        "could make estimate time longer than expected.",
        "\nDo you want to continue?",
        "\ny/n"
      ),
      "n"
    )

    .small_number_warning(
      data_$obs_size,
      if (risk_estimate_method == "cross_validation") 30 else 0,
      paste0(
        "The size of the observational dataset is insufficient ",
        "for cross vlidation to yield optimal results.",
        "\nSUGGESTION: Please consider utilizing alternative methods",
        "(e.g. risk_estimate_method = parametric_bootstrap) ",
        "for better estimate."
      )
    )

    risk_estimate_collect <- c()
    optimal_tau <- NULL

    for (tau in taus) {
      risk_estimate_collect <- append(
        risk_estimate_collect,
        do.call(
          risk_estimate_method,
          list(
            formula = formula,
            family = family,
            data_ = data_,
            tau = tau,
            discrepancy_function = discrepancy_function,
            tau_0 = tau_0,
            bootstrap_iteration_times = bootstrap_iteration_times,
            cv_fold_num = cv_fold_num
          )
        )
      )
    }

    optimal_tau <- taus[which.min(risk_estimate_collect)]

    glm_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = data_$data,
      weights = c(
        rep(1, data_$obs_size),
        rep(optimal_tau / data_$syn_size, data_$syn_size)
      )
    ))

    .large_number_warning(
      summary(glm_model)$coefficients[, grepl(
        "Pr.*",
        colnames(summary(glm_model)$coefficients)
      )],
      0.05,
      paste0(
        "Some(one) of the coefficient estimates in the model are(is) ",
        "not statistically significant at the ",
        "conventional significance level of 0.05.\n",
        "This suggests that there may not be sufficient evidence ",
        "to conclude that the variable(s) has a significant effect ",
        "on the response variables.\n",
        "SUGGESTION: Please consider utilizing althernative `risk_estimate_method` ",
        "or `catalytic::cat_glm_bayesian_joint` for better estimate."
      )
    )

    cat_glm_obj <- list(
      function_name = "cat_glm_tune",
      formula = formula,
      family = family,
      data_ = data_,
      tau = optimal_tau,
      coefficients = glm_model$coefficients,
      model = glm_model,
      min_risk_estimate = min(risk_estimate_collect),
      taus = taus,
      risk_estimate_collect = risk_estimate_collect,
      discrepancy_function = discrepancy_function,
      risk_estimate_method = risk_estimate_method,
      tau_0 = tau_0,
      bootstrap_iteration_times = bootstrap_iteration_times,
      cv_fold_num = cv_fold_num
    )

    class(cat_glm_obj) <- c(cat_glm_obj$class, "none_bayes", "cat_glm_tune")

    return(cat_glm_obj)
  }


parametric_bootstrap <- function(formula,
                                 family = gaussian(),
                                 data_,
                                 tau,
                                 ...,
                                 discrepancy_function = "square_error",
                                 tau_0 = 1,
                                 bootstrap_iteration_times = 100) {
  fam <- get_family(family)

  beta_tau_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau / data_$syn_size, data_$syn_size)
    )
  ))


  D <- mean(get_discrepancy(discrepancy_function)(
    X = data_$adj_obs_x,
    Y = data_$obs_y,
    model = beta_tau_model,
    fam = fam)[[1]])

  hat_tau_0_beta_model <- do.call(
    stats::glm,
    list(
      formula = formula,
      family = family,
      data = data_$data,
      weights = c(
        rep(1, data_$obs_size),
        rep(tau_0 / data_$syn_size, data_$syn_size)
      )
    )
  )

  mean_Y <- get_mean(fam)(
    X = data_$adj_obs_x,
    model = hat_tau_0_beta_model)

  sd_Y <- sd(mean_Y)

  boot_Y_collect <-
    lambda_boot_hat_mu_collect <- matrix(NA,
      nrow = data_$obs_size,
      ncol = bootstrap_iteration_times
    )

  for (bootstrap_iter in 1:bootstrap_iteration_times) {
    boot_Y <- get_data(
      fam
    )(n = data_$obs_size,
      mean = mean_Y,
      sd = sd_Y
    )
    boot_Y_df <- as.data.frame(boot_Y)
    names(boot_Y_df) <- names(data_$syn_y)

    boot_hat_tau_beta_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = cbind(
        data_$x,
        rbind(boot_Y_df, data_$syn_y)
      ),
      weights = c(
        rep(1, data_$obs_size),
        rep(tau / data_$syn_size, data_$syn_size)
      )
    ))

    boot_Y_collect[, bootstrap_iter] <- boot_Y

    lambda_boot_hat_mu_collect[, bootstrap_iter] <- lambda_function(
      discrepancy_function
    )(
      X = data_$adj_obs_x,
      coefficients = boot_hat_tau_beta_model$coefficients,
      fam = fam)
  }

  return(D + sum(sweep(boot_Y_collect, 1, rowMeans(boot_Y_collect)) *
    lambda_boot_hat_mu_collect) / bootstrap_iteration_times)
}

steinian_estimate <- function(formula,
                              family = gaussian(),
                              data_,
                              tau,
                              ...,
                              tau_0 = 1) {
  fam <- get_family(family)

  .not_same_stop(
    fam,
    "binomial",
    paste(
      "Please note that only logistic regression can use,",
      "`risk_estimate_method = steinian_estimate`."
    )
  )


  beta_tau_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau / data_$syn_size, data_$syn_size)
    )
  ))

  D <- mean(get_discrepancy(
    "logistic_deviance"
  )(
    X = data_$adj_obs_x,
    Y = data_$obs_y,
    model = beta_tau_model,
    fam)[[1]])

  hat_beta_0_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau_0 / data_$syn_size, data_$syn_size)
    )
  ))

  hat_mu_0 <- get_mean(fam)(X = data_$adj_obs_x,
    coefficients = hat_beta_0_model$coefficients)

  with_ <- with(hat_beta_0_model, fitted.values * (1 - fitted.values))

  steinian_sum <- 0
  for (obs in 1:data_$obs_size) {
    triangle_response <- data_$y
    triangle_response[obs, ] <- 1 - data_$y[obs, ]

    trigangle_hat_tau_beta_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = cbind(
        data_$adj_x,
        triangle_response
      ),
      weights = c(
        rep(1, data_$obs_size),
        rep(tau / data_$syn_size, data_$syn_size)
      )
    ))

    steinian_sum <- steinian_sum +
      with_[obs] * (2 * data_$y[obs, ] - 1) * (tcrossprod(
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
                             ...,
                             tau_0 = 1) {
  fam <- get_family(family)

  .not_same_stop(
    fam,
    "gaussian",
    paste(
      "Please note that only linear regression can use",
      "`risk_estimate_method = mallows_estimate`"
    )
  )

  beta_tau_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau / data_$syn_size, data_$syn_size)
    )
  ))

  D <- sum(get_discrepancy("square_error")(
    data_$adj_obs_x,
    data_$obs_y,
    beta_tau_model,
    fam)[[1]])

  beta_0_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau_0 / data_$syn_size, data_$syn_size)
    )
  ))


  var_Y <- sum((get_sim_y(fam)(X = data_$adj_obs_x,
    model = beta_0_model) -
    data_$obs_y)^2) / data_$obs_size

  tau_H <- H_matrix(
    X = data_$adj_obs_x,
    X_star = data_$adj_syn_x,
    tau = tau,
    M = data_$syn_size,
    n = data_$obs_size
  )

  return(D + 2 * var_Y * sum(diag(tau_H)) / data_$obs_size)
}

cross_validation <- function(formula,
                             family = gaussian(),
                             data_,
                             tau,
                             ...,
                             discrepancy_function = "square_error",
                             cv_fold_num = 5) {
  error_sum <- 0
  random_sample <- sample(1:data_$obs_size)
  fam <- get_family(family)

  for (fold in 1:cv_fold_num) {
    test_x <- data_$adj_obs_x[random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]
    test_y <- data_$obs_y[random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]

    train_x <- data_$adj_obs_x[-random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]
    y <- data_$obs_y[-random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)], ]

    estimate_beta_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = rbind(cbind(train_x, y), data_$syn_data),
      weights = c(
        rep(1, nrow(train_x)),
        rep(tau / data_$syn_size, data_$syn_size)
      )
    ))

    error_sum <- error_sum + mean(get_discrepancy(discrepancy_function)(
      X = test_x,
      Y = test_y,
      model = estimate_beta_model,
      fam = fam)[[1]])
  }

  return(error_sum / cv_fold_num)
}

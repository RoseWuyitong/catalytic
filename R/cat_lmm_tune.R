#' Title
#'
#' @param data_
#' @param taus
#' @param ...
#' @param discrepancy_function
#' @param cv_fold_num
#' @param risk_estimate_method
#'
#' @return
#' @export
#'
#' @examples
cat_lmm_tune <- function(data_,
                         taus = NULL,
                         ...,
                         discrepancy_function = "square_error",
                         cv_fold_num = 5,
                         risk_estimate_method = "cross_validation_lmm") {
  if (is.null(taus)) {
    taus <- seq(
      max(0, ncol(data_$obs_x) / 4 - 2),
      ncol(data_$obs_x) / 4 + 2, 0.5
    )
  }

  risk_estimate_collect <- do.call(
    risk_estimate_method,
    list(
      data_ = data_,
      tau = taus,
      discrepancy_function = discrepancy_function,
      cv_fold_num = cv_fold_num,
      ...
    )
  )

  optimal_tau <- taus[which.min(risk_estimate_collect)]
  optimal_model <- do.call(lme4::lmer, list(
    formula = data_$y ~ 1 + data_$x + (1 | data_$group_col),
    weight = c(
      rep(1, data_$obs_size),
      rep(optimal_tau / data_$syn_size, data_$syn_size)
    ),
    ...
  ))

  cat_lmm_tune_obj <- list(
    function_name = "cat_lmm_tune",
    data_ = data_,
    taus = taus,
    discrepancy_function = discrepancy_function,
    cv_fold_num = cv_fold_num,
    risk_estimate_method = risk_estimate_method,
    risk_estimate_collect = risk_estimate_collect,
    tau = optimal_tau,
    model = optimal_model,
    min_risk_estimate = min(risk_estimate_collect),
    ...
  )

  class(cat_lmm_tune_obj) <- c(cat_lmm_tune_obj$class, "cat_lmm_tune")

  return(cat_lmm_tune_obj)
}

cross_validation_lmm <- function(data_,
                                 taus,
                                 ...,
                                 discrepancy_function = "square_error",
                                 cv_fold_num = 5) {
  random_sample <- sample(1:data_$obs_size)

  risk_estimate_collect <- rep(0, length(taus))
  for (fold in 1:cv_fold_num) {
    idx <- random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)]

    test_x <- data_$obs_x[idx, ]
    test_y <- as.matrix(data_$obs_y[idx, ])
    test_z <- data_$obs_z[idx, ]
    train_group <- as.matrix(data_$obs_group[idx, ])

    train_x <- data_$obs_x[-idx, ]
    train_y <- as.matrix(data_$obs_y[-idx, ])
    train_z <- data_$obs_z[-idx, ]
    train_group <- as.matrix(data_$obs_group[-idx, ])

    risk_estimate_collect_tau <- c()
    for (tau in taus) {
      yy <- rbind(train_y, data_$syn_y)
      xx <- rbind(train_x, data_$syn_x)
      gg <- rbind(train_group, data_$syn_group_col)
      cv_model <- do.call(lme4::lmer, list(
        formula = yy ~ xx + (1 | gg),
        weight = c(rep(1, nrow(train_x)), rep(tau / data_$syn_size, data_$syn_size))
      ))

      error_sum <- mean(get_discrepancy(discrepancy_function)(
        X = test_x,
        Y = test_y,
        est_coefficients = lme4::fixef(cv_model),
        fam = "gaussian")[[1]])

      risk_estimate_collect_tau <- c(risk_estimate_collect_tau, error_sum / cv_fold_num)
    }
    risk_estimate_collect <- risk_estimate_collect + risk_estimate_collect_tau
  }

  return(risk_estimate_collect)
}

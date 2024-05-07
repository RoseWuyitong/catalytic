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

  if (is.null(taus)) taus <- seq(max(0, ncol(data_$obs_x) / 4 - 2), ncol(data_$data) / 4 + 2, 0.5)

  risk_estimate_collect <- matrix(nrow = 0,
                                  ncol = 2,
                                  dimnames = list(NULL, c("our_error", "lmer_error")))

  for (tau in taus) {
    risk_estimate_collect <- rbind(
      risk_estimate_collect,
      do.call(
        risk_estimate_method,
        list(
          data_ = data_,
          tau = tau,
          discrepancy_function = discrepancy_function,
          cv_fold_num = cv_fold_num,
          ...
        )
      )
    )
  }

  our_optimal_tau <- taus[which.min(risk_estimate_collect[,"our_error"])]
  our_optimal_model <- do.call(cat_lmm, list(
    data_ = data_,
    tau = our_optimal_tau,
    ...
  ))

  lmer_optimal_tau <- taus[which.min(risk_estimate_collect[,"lmer_error"])]
  lmer_optimal_model <- do.call(lme4::lmer, list(
    formula=data_$y ~ 1 + data_$x + ( 1 | data_$group_col),
    weight=c(rep(1, data_$obs_size), rep(lmer_optimal_tau/data_$syn_size, data_$syn_size)
    )))


  cat_lmm_tune_obj <- list(
    function_name = "cat_lmm_tune",
    data_ = data_,
    taus = taus,
    discrepancy_function = discrepancy_function,
    cv_fold_num = cv_fold_num,
    risk_estimate_method = risk_estimate_method,
    risk_estimate_collect = risk_estimate_collect,
    ...,
    our_tau = our_optimal_tau,
    our_model = our_optimal_model,
    our_min_risk_estimate = min(risk_estimate_collect[,"our_error"]),
    lmer_tau = lmer_optimal_tau,
    lmer_model = lmer_optimal_model,
    lmer_min_risk_estimate = min(risk_estimate_collect[,"lmer_error"])
  )

  class(cat_lmm_tune_obj) <- c(cat_lmm_tune_obj$class, "cat_lmm_tune")

  return(cat_lmm_tune_obj)
}

cross_validation_lmm <- function(data_,
                                 tau,
                                 ...,
                                 discrepancy_function = "square_error",
                                 cv_fold_num = 5) {
  error_sum <- lmer_error_sum <- 0
  random_sample <- sample(1:data_$obs_size)

  for (fold in 1:cv_fold_num) {
    idx <- random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
                                (data_$obs_size * fold / cv_fold_num)]

    test_x <- data_$obs_x[idx, ]
    test_y <- data_$obs_y[idx, ]
    test_z <- data_$obs_z[idx, ]
    train_x <- data_$obs_x[-idx, ]
    train_y <- data_$obs_y[-idx, ]
    train_z <- data_$obs_z[-idx, ]

    cv_data_ <- data_
    cv_data_$obs_x <- train_x
    cv_data_$obs_y <- train_y
    cv_data_$obs_z <- train_z

    estimate_model <- do.call(cat_lmm, list(
      data_ = cv_data_,
      tau = tau,
      ...
    ))

    lmer_model <- do.call(lme4::lmer, list(
      formula=data_$y ~ 1 + data_$x + ( 1 | data_$group_col),
      weight=c(rep(1, data_$obs_size), rep(tau/data_$syn_size, data_$syn_size)
    )))

    error_sum <- error_sum + mean(get_discrepancy(discrepancy_function)(
      X = test_x,
      Y = test_y,
      est_coefficients = estimate_model$coefficients,
      fam = "gaussian")[[1]])

    lmer_error_sum <- lmer_error_sum + mean(get_discrepancy(discrepancy_function)(
      X = test_x,
      Y = test_y,
      est_coefficients = lme4::fixef(lmer_model),
      fam = "gaussian")[[1]])
  }

  return(c(error_sum / cv_fold_num,lmer_error_sum / cv_fold_num))
}

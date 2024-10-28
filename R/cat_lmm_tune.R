#' Title
#'
#' @param cat_init
#' @param tau_seq
#' @param ...
#' @param "square_error"
#' @param cross_validation_fold_num
#' @param "cross_validation_lmm"
#'
#' @return
#' @export
#'
#' @examples
cat_lmm_tune <- function(
    cat_init,
    tau_seq = NULL,
    cross_validation_fold_num = 5) {
  # Validate Input Parameters
  validate_lmm_input(
    cat_init = cat_init,
    tau_seq = tau_seq,
    cross_validation_fold_num = cross_validation_fold_num,
  )

  if (is.null(tau_seq)) {
    tau_seq <- seq(
      max(0, ncol(cat_init$adj_x) / 4 - 2),
      ncol(cat_init$adj_x) / 4 + 2, 0.5
    )
  }

  risk_estimate_list <- do.call(
    "cross_validation_lmm",
    list(
      cat_init = cat_init,
      tau_seq = tau_seq,
      cross_validation_fold_num = cross_validation_fold_num
    )
  )

  optimal_tau <- tau_seq[which.min(risk_estimate_list)]

  optimal_model <- do.call(lme4::lmer, list(
    formula = cat_init$y ~ 1 + cat_init$x + (1 | cat_init$group_col),
    weight = c(
      rep(1, cat_init$obs_size),
      rep(optimal_tau / cat_init$syn_size, cat_init$syn_size)
    )
  ))

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_lmm_tune",
    ## Input/Processed parameters
    cat_init = cat_init,
    tau_seq = tau_seq,
    cross_validation_fold_num = cross_validation_fold_num,
    ## Result
    tau = optimal_tau,
    model = optimal_model,
    coefficients = coef(optimal_model),
    risk_estimate_list = risk_estimate_list
  )

  class(cat_model) <- c(cat_model$class, "cat_lmm_tune", "cat_tune", "cat_lmm")

  return(cat_model)
}

cross_validation_lmm <- function(
    cat_init,
    tau_seq,
    cross_validation_fold_num = 5) {
  random_sample <- sample(1:cat_init$obs_size)
  risk_estimate_list <- rep(0, length(tau_seq))

  # The perform of cross-validation
  for (fold in seq_len(cross_validation_fold_num)) {
    test_idx <- random_sample[
      (cat_init$obs_size * (fold - 1) / cross_validation_fold_num + 1):
      (cat_init$obs_size * fold / cross_validation_fold_num)
    ]

    test_x <- cat_init$obs_x[test_idx, , drop = FALSE]
    test_y <- cat_init$obs_y[test_idx]
    test_z <- cat_init$obs_z[test_idx, , drop = FALSE]
    test_group <- cat_init$obs_group[test_idx]

    train_x <- cat_init$obs_x[-test_idx, , drop = FALSE]
    train_y <- cat_init$obs_y[-test_idx]
    train_z <- cat_init$obs_z[-test_idx, , drop = FALSE]
    train_group <- cat_init$obs_group[-test_idx]

    risk_estimate_list_fold <- c()

    # The iteration over the tau sequence for estimation
    for (tau in tau_seq) {
      cv_y <- rbind(train_y, cat_init$syn_y)
      cv_x <- rbind(train_x, cat_init$syn_x)
      cv_group <- rbind(train_group, cat_init$syn_group_col)

      cv_model <- do.call(
        lme4::lmer, list(
          formula = cv_y ~ cv_x + (1 | cv_group),
          weight = c(
            rep(1, nrow(train_x)),
            rep(tau / cat_init$syn_size, cat_init$syn_size)
          )
        )
      )
      browser()
      # TODO: check use predict(model, newData)
      prediction_error <- mean(
        get_discrepancy("square_error")(
          family_string = "gaussian",
          X = test_x,
          Y = test_y,
          coefs = lme4::fixef(cv_model))
      )

      risk_estimate_list_fold <- c(
        risk_estimate_list_fold,
        prediction_error
      )
    }
    risk_estimate_list <- risk_estimate_list + risk_estimate_list_fold
  }

  return(risk_estimate_list / cross_validation_fold_num)
}

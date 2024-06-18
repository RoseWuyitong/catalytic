#' Title
#'
#' @param data_
#' @param taus
#' @param ...
#' @param discrepancy_function
#' @param cv_fold_num
#' @param likelihood_method
#'
#' @return
#' @export
#'
#' @examples
cat_cox_tune <- function(
    formula,
    data_,
    taus = NULL,
    method = c("CRE", "WME"),
    ...,
    cv_fold_num = 5,
    likelihood_method = "cross_validation_cox") {
  .not_true_stop(
    is.data_cox(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization_cox."
    )
  )

  # Extract time, event and other columns(x) from `data`
  f_lhs <- get_lhs(formula)

  if (!is.null(f_lhs)) {
    .not_true_stop(
      is.surv_formula(formula),
      "Please use `Surv` for formula."
    )
    time_col_name <- gsub("`", "", f_lhs[2])
    event_col_name <- gsub("`", "", f_lhs[3])

    .not_true_stop(
      time_col_name == data_$time_col_name & event_col_name == data_$event_col_name,
      "Please use same time and event column name for formula"
    )
  }

  if (is.null(taus)) {
    taus <- seq(
      max(0, ncol(data_$obs_x) / 4 - 2),
      ncol(data_$obs_x) / 4 + 2, 0.5
    )
  }

  method <- match.arg(method)

  likelihood_collect <- do.call(
    likelihood_method,
    list(
      formula = formula,
      data_ = data_,
      method = method,
      taus = taus,
      cv_fold_num = cv_fold_num,
      entry_points = data_$entry_points,
      ...
    )
  )

  optimal_tau <- taus[which.max(likelihood_collect)]

  optimal_model <- do.call(cat_cox, list(
    formula = formula,
    data_ = data_,
    method = method,
    tau = optimal_tau
  ))

  cat_cox_tune_obj <- list(
    function_name = "cat_cox_tune",
    method = method,
    data_ = data_,
    taus = taus,
    cv_fold_num = cv_fold_num,
    likelihood_method = likelihood_method,
    likelihood_collect = likelihood_collect,
    tau = optimal_tau,
    model = optimal_model,
    coefficients = coef(optimal_model),
    max_likelihood = max(likelihood_collect),
    ...
  )

  class(cat_cox_tune_obj) <- c(cat_cox_tune_obj$class, "cat_cox_tune")

  return(cat_cox_tune_obj)
}

cross_validation_cox <- function(formula,
                                 data_,
                                 method,
                                 taus,
                                 ...,
                                 discrepancy_function = "square_error",
                                 cv_fold_num = 5) {
  random_sample <- sample(1:data_$obs_size)

  likelihood_collect <- rep(0, length(taus))

  for (fold in 1:cv_fold_num) {
    idx <- random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):(data_$obs_size * fold / cv_fold_num)]

    train_x <- data_$obs_x[-idx, ]
    train_time <- data_$obs_time[-idx, ]
    train_event <- data_$obs_event[-idx, ]

    cv_data_ <- data_
    cv_data_$obs_size <- nrow(train_x)
    cv_data_$obs_x <- as.matrix(train_x)
    cv_data_$obs_time <- matrix(train_time,
      ncol = 1,
      dimnames = list(NULL, data_$time_col_name)
    )
    cv_data_$obs_event <- matrix(train_event,
      ncol = 1,
      dimnames = list(NULL, data_$event_col_name)
    )

    cv_data_$obs_data <- cbind(
      cv_data_$obs_x,
      cv_data_$obs_time,
      cv_data_$obs_event
    )
    cv_data_$x <- rbind(cv_data_$obs_x, cv_data_$syn_x)
    cv_data_$time <- rbind(cv_data_$obs_time, cv_data_$syn_time)
    cv_data_$event <- rbind(cv_data_$obs_event, cv_data_$syn_event)
    cv_data_$size <- cv_data_$obs_size + cv_data_$syn_size
    cv_data_$entry_points <- cv_data_$entry_points[-idx]
    likelihood_collect_tau <- c()

    for (tau in taus) {
      cv_model <- do.call(
        cat_cox,
        list(
          formula = formula,
          data_ = cv_data_,
          method = method,
          tau = tau
        )
      )

      error_sum <- get_cox_partial_likelihood(
        X = data_$obs_x,
        time = data_$obs_time,
        event = data_$obs_event,
        coefs = cv_model$coefficients,
        entry_points = data_$entry_points
      ) - get_cox_partial_likelihood(
        X = train_x,
        time = train_time,
        event = train_event,
        coefs = cv_model$coefficients,
        entry_points = cv_data_$entry_points
      )

      likelihood_collect_tau <- c(
        likelihood_collect_tau,
        error_sum / cv_fold_num
      )
    }

    likelihood_collect <- likelihood_collect + likelihood_collect_tau
  }

  return(likelihood_collect)
}

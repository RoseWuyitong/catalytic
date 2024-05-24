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
cat_cox_tune <- function(data_,
                         taus = NULL,
                         method = c("CRE", "WME"),
                         ...,
                         cv_fold_num = 5,
                         likelihood_method = "cross_validation_cox",
                         hazard_value = NULL) {
  if (is.null(taus)) {
    taus <- seq(
      max(0, ncol(data_$obs_x) / 4 - 2),
      ncol(data_$obs_x) / 4 + 2, 0.5
    )
  }
  if (is.null(hazard_value)) hazard_value <- sum(data_$obs_event) / mean(data_$obs_time) / data_$obs_size
  method <- match.arg(method)

  likelihood_collect <- do.call(
    likelihood_method,
    list(
      data_ = data_,
      method = method,
      taus = taus,
      cv_fold_num = cv_fold_num,
      ...
    )
  )

  optimal_tau <- taus[which.max(likelihood_collect)]

  optimal_model <- do.call(cat_cox, list(
    data_ = data_,
    method = method,
    tau = optimal_tau,
    ...
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

cross_validation_cox <- function(data_,
                                 method,
                                 taus,
                                 ...,
                                 discrepancy_function = "square_error",
                                 cv_fold_num = 5) {
  # TODO: user defined
  entry <- rep(0, data_$obs_size)

  random_sample <- sample(1:data_$obs_size)

  likelihood_collect <- rep(0, length(taus))

  for (fold in 1:cv_fold_num) {
    idx <- random_sample[(data_$obs_size * (fold - 1) / cv_fold_num + 1):
    (data_$obs_size * fold / cv_fold_num)]

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
    )[, colnames(data_$obs_data)]
    cv_data_$x <- rbind(cv_data_$obs_x, cv_data_$syn_x)
    cv_data_$time <- rbind(cv_data_$obs_time, cv_data_$syn_time)
    cv_data_$event <- rbind(cv_data_$obs_event, cv_data_$syn_event)

    likelihood_collect_tau <- c()

    for (tau in taus) {
      cat("\n fold: ", fold, "cv tau:", tau)
      cv_model <- do.call(
        cat_cox,
        list(
          data_ = cv_data_,
          method = method,
          tau = tau,
          ...
        )
      )

      error_sum <- cox_pl(
        cv_model$coefficients,
        data_$obs_x,
        rep(0, nrow(test_x)),
        data_$obs_time,
        data_$obs_event
      ) - cox_pl(
        cv_model$coefficients,
        train_x,
        rep(0, nrow(train_x)),
        train_time,
        train_event
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

cox_pl <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  entry <- rep(0, length(Ts))
  uncersor_set <- which(as.numeric(event) == 1)
  each_PL <- sapply(uncersor_set, function(i) {
    ts <- Ts[i]
    xs <- Xs[i, ]
    riskindex <- get_risk_set(ts, entry, Ts, event)
    if (length(riskindex) > 1) {
      X_risk_set <- Xs[riskindex, ]
    } else {
      X_risk_set <- matrix(Xs[riskindex, ], nc = ncol(Xs))
    }
    # X_risk_set <- Xs[riskindex,] %>%
    #  matrix(ncol = ncol(Xs))
    theta <- exp(X_risk_set %*% beta)
    return(sum(xs * beta) - log(sum(theta)))
  })

  return(sum(each_PL))
}

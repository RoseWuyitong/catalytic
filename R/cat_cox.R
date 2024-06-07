#' Title
#'
#' @param formula
#' @param data_
#' @param tau
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cat_cox <-
  function(formula,
           data_,
           tau = NULL,
           method = c("CRE", "WME"),
           ...,
           initial_values = NULL,
           tol = 1e-5,
           max_iter = 25) {
    .not_true_stop(
      is.data_cox(data_),
      paste(
        "Please enter `data_` that generated from",
        "catalytic::data_initialization_cox."
      )
    )

    if (is.null(tau)) tau <- ncol(data_$obs_x) / 4

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

    method <- match.arg(method)

    data_ <- get_adj_data(data_ = data_, formula = formula)

    if (method == "CRE") {
      coefs <- if (is.null(initial_values)) rep(0, ncol(data_$adj_obs_x)) else initial_values


      diff_new_old <- TRUE
      iter <- 0
      iter_collect <- matrix(
        data = c(coefs),
        ncol = (length(coefs)),
        dimnames = list(NULL, paste0("beta[", 1:length(coefs), "]"))
      )

      coefs_old <- coefs
      coef_matrix <- get_cox_hessian(
        X = data_$adj_obs_x,
        time = data_$obs_time,
        event = data_$obs_event,
        coefs = coefs_old,
        entry_points = data_$entry_points
      ) + tau / data_$syn_size * get_syn_cox_hessian(
        X = data_$adj_syn_x,
        time = data_$syn_time,
        coefs = coefs_old,
        # HERE
        hazard_constant = data_$hazard_constant
      )

      constant_vector <- get_cox_gradient(
        X = data_$adj_obs_x,
        time = data_$obs_time,
        event = data_$obs_event,
        coefs = coefs_old,
        entry_points = data_$entry_points
      ) + tau / data_$syn_size * get_syn_cox_gradient(
        X = data_$adj_syn_x,
        time = as.vector(data_$syn_time),
        coefs = coefs_old,
        hazard_constant = data_$hazard_constant
      )

      coefs <- coefs_old - get_qr_solve(coef_matrix, constant_vector)

      while ((iter < max_iter) & diff_new_old) {
        coefs_old <- coefs
        coef_matrix <- get_cox_hessian(
          X = data_$adj_obs_x,
          time = data_$obs_time,
          event = data_$obs_event,
          coefs = coefs_old,
          entry_points = data_$entry_points
        ) + tau / data_$syn_size * get_syn_cox_hessian(
          X = data_$adj_syn_x,
          time = data_$syn_time,
          coefs = coefs_old,
          hazard_constant = sum(data_$syn_event) / mean(data_$syn_time) / nrow(data_$syn_x)
        )

        constant_vector <- get_cox_gradient(
          X = data_$adj_obs_x,
          time = data_$obs_time,
          event = data_$obs_event,
          coefs = coefs_old,
          entry_points = data_$entry_points
        ) + tau / data_$syn_size * get_syn_cox_gradient(
          X = data_$adj_syn_x,
          time = as.vector(data_$syn_time),
          coefs = coefs_old,
          hazard_constant = sum(data_$syn_event) / mean(data_$syn_time) / nrow(data_$syn_x)
        )

        coefs <- coefs_old - get_qr_solve(coef_matrix, constant_vector)

        iter <- iter + 1
        iter_collect <- rbind(iter_collect, c(coefs))

        .true_warning(
          any(is.na(coefs)),
          "Please note NA appears in iteration."
        )

        .true_warning(
          norm(as.numeric(coefs), "2") > 1000000,
          "Please not the euclidean norm of coefficients is too large."
        )

        diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > tol * norm(as.numeric(coefs_old), "2")
      }

      .true_warning(
        iter == max_iter,
        "Please not it reached maximun iteration, it is not converging."
      )
    } else {
      iter_collect <- NULL
      iter <- NULL
      coefs_model <- do.call(survival::coxph, list(
        formula = survival::Surv(data_$time, data_$event) ~ data_$adj_x,
        weights = c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size)),
        ties = "breslow"
      ))

      coefs <- coef(coefs_model)
    }


    data_ <- list(
      function_name = "cat_cox",
      method = method,
      data_ = data_,
      tau = tau,
      iter = iter,
      iter_collect = iter_collect,
      coefficients = as.numeric(coefs),
      initial_values = initial_values,
      hazard_constant = hazard_constant,
      tol = tol,
      max_iter = max_iter,
      ...
    )

    class(data_) <- c(data_$class, "cat_cox")

    return(data_)
  }

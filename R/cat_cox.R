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
  function(data_,
           tau = NULL,
           method = c("CRE", "WME"),
           ...,
           initial_values = NULL,
           hazard_function = NULL,
           tol = 1e-5,
           max_iter = 25,
           entry = NULL) {
    .not_true_stop(
      is.data_cox(data_),
      paste(
        "Please enter `data_` that generated from",
        "catalytic::data_initialization_cox."
      )
    )

    if (is.null(tau)) tau <- ncol(data_$obs_x) / 4
    method <- match.arg(method)
    if (method == "CRE") {
      if (is.null(hazard_function)) hazard_function <- function(t) sum(data_$obs_event) / mean(data_$obs_time) / data_$obs_size
      if (is.null(entry)) entry <- rep(0, data_$obs_size)
      coefs <- if (is.null(initial_values)) rep(0, ncol(data_$obs_x)) else initial_values


      diff_new_old <- TRUE
      iter <- 0
      iter_collect <- matrix(
        data = c(coefs),
        ncol = (length(coefs)),
        dimnames = list(NULL, paste0("beta[", 1:length(coefs), "]"))
      )

      while ((iter < max_iter) & diff_new_old) {
        coefs_old <- coefs
        syn_ln <- as.numeric(data_$syn_x %*% coefs_old)

        coef_matrix <- get_cox_hessian(
          coefs_old,
          data_$obs_x,
          entry,
          data_$obs_time,
          data_$obs_event
        ) + tau / data_$syn_size * get_syn_cox_hessian(
          syn_ln,
          data_$syn_x,
          data_$syn_time,
          hazard_function
        )

        constant_vector <- get_cox_gradient(
          coefs,
          data_$obs_x,
          entry,
          data_$obs_time,
          data_$obs_event
        ) + tau / data_$syn_size * get_syn_cox_gradient(
          syn_ln,
          data_$syn_x,
          data_$syn_time,
          data_$syn_event,
          hazard_function
        )

        coefs <- coefs_old - get_qr_solve(coef_matrix, constant_vector)

        iter <- iter + 1
        iter_collect <- rbind(iter_collect, c(coefs))

        if (is.na(coefs[1])) {
          warning("NA BLEAK UP")
          return(initial_values + 10000)
        }
        if (norm(as.numeric(coefs), "2") > 1000000) {
          warning("BLEAK UP")
          return(initial_values + 50000)
        }
        if (norm(as.numeric(coefs_old - coefs), "2") > tol * norm(as.numeric(coefs_old), "2")) {
          diff_new_old <- TRUE
        } else {
          diff_new_old <- FALSE
        }
      }

      if (iter == max_iter) {
        warning("Reached maximun iteration, please note it is not converging.")
      }
    } else {
      iter_collect <- NULL
      iter <- NULL

      coefs_model <- do.call(survival::coxph, list(
        formula = Surv(data_$time, data_$event) ~ data_$x,
        weights = c(rep(1, data_$obs_size), rep(tau / data_$syn_size, data_$syn_size)),
        ties = "breslow",
        ...
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
      hazard_function = hazard_function,
      tol = tol,
      max_iter = max_iter,
      entry = entry,
      ...
    )
    class(data_) <- c(data_$class, "cat_cox")

    return(data_)
  }

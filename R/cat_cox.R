#' Fit Cox Proportional-Hazards Model with Synthetic Data
#'
#' This function fits a Cox Proportional-Hazards Model using the specified formula and data.
#' The function supports two methods for fitting the model: "CRE" (Cox Regression Estimator) and "WME" (Weighted Mean Estimator).
#'
#' @param formula The formula for the Cox model, specifying the time and status variables and covariates.
#' @param cat_init The initialized data object created from `catalytic::data_initialization_cox`.
#' @param tau A regularization parameter. Default is `NULL`, which sets tau to the number of covariates.
#' @param method The method for fitting the Cox model, either "CRE" (Cox Regression Estimator) or "WME" (Weighted Mean Estimator). Default is "CRE".
#' @param ... Additional arguments passed to the Cox model function.
#' @param init_coefficients Initial coefficients for the Cox model. Default is `NULL`.
#' @param tol The tolerance for convergence. Default is 1e-5.
#' @param max_iter The maximum number of iterations for convergence. Default is 25.
#'
#' @return A list containing the fitted Cox model and various fitting details:
#' \describe{
#'   \item{function_name}{The name of the function, `cat_cox`.}
#'   \item{formula}{The formula used for the Cox model.}
#'   \item{cat_init}{The initialized data object used for the model.}
#'   \item{tau}{The regularization parameter used.}
#'   \item{method}{The method used for fitting the model, either "CRE" or "WME".}
#'   \item{coefficients}{The coefficients of the fitted Cox model.}
#'   \item{model}{The fitted survival model (only for "WME" method).}
#'   \item{iter}{The number of iterations performed (only for "CRE" method).}
#'   \item{iteration_log}{The collected coefficients for each iteration (only for "CRE" method).}
#'   \item{init_coefficients}{The initial coefficients used (if any).}
#'   \item{tol}{The tolerance for convergence.}
#'   \item{max_iter}{The maximum number of iterations allowed.}
#' }
#' @export
#'
#' @examples
#' # Assuming `data_cox` is an initialized data object created with `data_initialization_cox`
#' formula <- survival::Surv(time, status) ~ .
#' cat_cox_model <- cat_cox(
#'   formula = formula,
#'   data_ = data_cox,
#'   tau = 5,
#'   method = "CRE",
#'   max_iter = 25
#' )
#'
#' print(cat_cox_model)
cat_cox <- function(formula,
                    cat_init,
                    tau = NULL,
                    method = c("CRE", "WME"),
                    init_coefficients = NULL,
                    tol = 1e-5,
                    max_iter = 25) {
  method <- match.arg(method)

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- as.formula(paste0(
    "survival::Surv(",
    cat_init$time_col_name,
    ",",
    cat_init$status_col_name,
    ")",
    f_rhs
  ))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  # Validate Input Parameters
  validate_cox_input(
    formula = formula,
    cat_init = cat_init,
    tau = tau,
    init_coefficients = init_coefficients,
    tol = tol,
    max_iter = max_iter
  )

  if (is.null(tau)) tau <- ncol(cat_init$adj_x)

  coefs_model <- NULL
  iteration_log <- NULL
  iter <- NULL

  if (method == "CRE") {
    coefs <- if (is.null(init_coefficients)) rep(0, ncol(cat_init$adj_obs_x)) else init_coefficients

    diff_new_old <- TRUE

    iter <- 0
    iteration_log <- matrix(
      data = c(coefs),
      ncol = (length(coefs)),
      dimnames = list(NULL, colnames(cat_init$adj_x))
    )

    while ((iter < max_iter) & diff_new_old) {
      coefs_old <- coefs
      hessian_matrix <- get_cox_hessian(
        X = cat_init$adj_obs_x,
        time = cat_init$obs_time,
        status = cat_init$obs_status,
        coefs = coefs_old,
        entry_points = cat_init$entry_points
      ) + tau / cat_init$syn_size * get_cox_syn_hessian(
        X = cat_init$adj_syn_x,
        time = cat_init$syn_time,
        coefs = coefs_old,
        hazard_constant = cat_init$hazard_constant
      )

      gradient_vector <- get_cox_gradient(
        X = cat_init$adj_obs_x,
        time = cat_init$obs_time,
        status = cat_init$obs_status,
        coefs = coefs_old,
        entry_points = cat_init$entry_points
      ) + tau / cat_init$syn_size * get_cox_syn_gradient(
        X = cat_init$adj_syn_x,
        time = cat_init$syn_time,
        coefs = coefs_old,
        hazard_constant = cat_init$hazard_constant
      )

      coefs <- coefs_old - get_cox_qr_solve(hessian_matrix, gradient_vector)

      iter <- iter + 1
      iteration_log <- rbind(iteration_log, c(coefs))

      diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > (tol * norm(as.numeric(coefs_old), "2"))
    }
  } else {
    coefs_model <- if (tau == 0) {
      do.call(survival::coxph, list(
        formula = full_formula,
        data = cat_init$obs_data,
        ties = "breslow"
      ))
    } else {
      do.call(survival::coxph, list(
        formula = full_formula,
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau / cat_init$syn_size, cat_init$syn_size)
        ),
        ties = "breslow"
      ))
    }

    coefs <- coef(coefs_model)
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_cox",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau = tau,
    method = method,
    init_coefficients = init_coefficients,
    tol = tol,
    max_iter = max_iter,
    ## Result
    iter = iter,
    iteration_log = iteration_log,
    model = coefs_model,
    coefficients = coefs
  )

  class(cat_model) <- c(cat_model$class, "cat_cox")

  return(cat_model)
}

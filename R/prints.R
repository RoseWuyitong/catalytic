#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cat_tau_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\ntau: ", x$tau,
    "\n\nEstimated beta: ", x$est_beta
  )
  invisible(x)
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
  cat("\nFormula for frequntist prior:  ")
  print(x$formula)
  cat(
    "\nRisk estimate method: ", x$risk_estimate_method,
    "\n\nOptimal tau: ", x$optimal_tau,
    "\n\nOptimal estimated beta: ", x$optimal_beta
  )
  invisible(x)
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
print.cat_bayesian_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nMCMC algorithm: ", x$mcmc_algorithm,
    "\n\nJoint tau: ", x$joint,
    "\n\nMean of estimated argments: ", x$arg_mean
  )
  invisible(x)
}

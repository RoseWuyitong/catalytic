#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.initialization <- function(x, ...) {
  cat("\nFormula for generating synthetic response:\n")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\nSynthetic data size: ", x$syn_size,
    "\n\nSynthetic data generation steps:\n\n"
  )
  print(x$syn_x_generation_steps)
  cat("\n")
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
print.cat_tau_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nFamily: ", x$family$family,
    "\n\ntau: ", x$tau,
    "\n\nEstimated beta: ", x$beta
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
    "\n\nOptimal tau: ", x$tau,
    "\n\nOptimal estimated beta: ", round(x$beta, 2)
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
  cat(
    "\nMCMC algorithm: ", x$mcmc_algorithm,
    "\n\nJoint tau: ", x$joint,
    "\n\nMean of estimated argments: ", x$beta
  )
  invisible(x)
}

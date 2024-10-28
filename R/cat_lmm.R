#' Title
#'
#' @param cat_init
#' @param tau
#' @param ...
#' @param max_iter
#' @param tol
#' @param variance_params_0
#' @param coefs_0
#'
#' @return
#' @export
#'
#' @examples
cat_lmm <- function(
    cat_init,
    tau = NULL,
    residual_variance_0 = 1,
    random_effect_variance_0 = 1,
    coefs_0 = NULL,
    optimize_domain = c(0, 30),
    max_iter = 500,
    tol = 1e-08) {
  # Validate Input Parameters
  validate_lmm_input(
    cat_init = cat_init,
    tau = tau,
    residual_variance_0 = residual_variance_0,
    random_effect_variance_0 = random_effect_variance_0,
    coefs_0 = coefs_0,
    optimize_domain = optimize_domain,
    max_iter = max_iter,
    tol = tol
  )

  if (is.null(tau)) tau <- ncol(cat_init$obs_x) / 4

  # Sort synthetic data by group
  sorted_syn_data <- cbind(
    cat_init$syn_x,
    setNames(as.data.frame(cat_init$syn_y), cat_init$y_col),
    cat_init$syn_z,
    setNames(as.data.frame(cat_init$syn_group), cat_init$group_col)
  )[order(cat_init$syn_group), ]

  cat_init$syn_data <- sorted_syn_data
  cat_init$syn_x <- sorted_syn_data$syn_x
  cat_init$syn_y <- sorted_syn_data$syn_y
  cat_init$syn_z <- sorted_syn_data$syn_z
  cat_init$syn_group <- sorted_syn_data$syn_group

  # Calculates the eigenvalues and eigenvectors of obs_z
  obs_z_eigen <- eigen(tcrossprod(as.matrix(cat_init$obs_z)))
  obs_z_eigenvalues <- obs_z_eigen$values
  obs_z_eigenvectors <- obs_z_eigen$vectors

  # Calculates the eigenvalues and eigenvectors of syb_z
  syn_z_eigenvalues <- c(rep(0, cat_init$syn_size))
  syn_z_eigenvectors <- matrix(0, ncol = cat_init$syn_size, nrow = cat_init$syn_size)
  row_offset <- 0

  for (g in unique(cat_init$group)) {
    group_idx <- which(cat_init$syn_group == g)

    syn_z_eigen <- eigen(tcrossprod(cat_init$syn_z[group_idx, ]))
    syn_z_eigenvalues[group_idx] <- syn_z_eigen$values
    syn_z_eigenvectors[((1:nrow(syn_z_eigen$vectors)) + row_offset), group_idx] <- syn_z_eigen$vectors

    row_offset <- row_offset + nrow(syn_z_eigen$vectors)
  }

  coefs <- if (is.null(coefs_0)) rnorm(ncol(cat_init$obs_x)) else coefs_0
  residual_variance <- residual_variance_0
  random_effect_variance <- random_effect_variance_0

  iteration_log <- matrix(
    data = c(residual_variance, random_effect_variance, coefs),
    ncol = (2 + length(coefs))
  )

  colnames(iteration_log) <- c(
    "residual_variance",
    "random_effect_variance",
    paste0("prev_coefs[", 1:length(coefs), "]")
  )

  iter <- 0
  parameters_converged <- FALSE

  while ((iter < max_iter) & (!parameters_converged)) {
    prev_coefs <- coefs
    prev_residual_variance <- residual_variance
    prev_random_effect_variance <- random_effect_variance

    obs_z_inv <- quadform::quad.tform(
      diag(1 / (residual_variance + random_effect_variance * obs_z_eigenvalues)),
      obs_z_eigenvectors
    )

    syn_z_inv <- quadform::quad.tform(
      diag(1 / (residual_variance + random_effect_variance * syn_z_eigenvalues)),
      syn_z_eigenvectors
    )

    # Update coefficients using generalized inverse and current parameter estimates
    coefs <- MASS::ginv(
      (tau / cat_init$syn_size) * quadform::quad.form(syn_z_inv, cat_init$syn_x) +
        quadform::quad.form(obs_z_inv, cat_init$obs_x)
    ) %*% (
      quadform::quad.3form(obs_z_inv, cat_init$obs_x, cat_init$obs_y) +
        tau / cat_init$syn_size * quadform::quad.3form(syn_z_inv, cat_init$syn_x, cat_init$syn_y))

    obs_adjusted_residuals <- t(obs_z_eigenvectors) %*% (cat_init$obs_y - cat_init$obs_x %*% coefs)
    syn_adjusted_residuals <- t(syn_z_eigenvectors) %*% (cat_init$syn_y - cat_init$syn_x %*% coefs)

    random_effect_variance <- stats::optimize(update_lmm_variance,
      optimize_domain,
      residual_variance = prev_residual_variance,
      obs_z_eigenvalues = obs_z_eigenvalues,
      syn_z_eigenvalues = syn_z_eigenvalues,
      obs_adjusted_residuals = obs_adjusted_residuals,
      syn_adjusted_residuals = syn_adjusted_residuals,
      tau = tau
    )$minimum

    residual_variance <- stats::optimize(update_lmm_variance,
      optimize_domain,
      random_effect_variance = random_effect_variance,
      obs_z_eigenvalues = obs_z_eigenvalues,
      syn_z_eigenvalues = syn_z_eigenvalues,
      obs_adjusted_residuals = obs_adjusted_residuals,
      syn_adjusted_residuals = syn_adjusted_residuals,
      tau = tau
    )$minimum

    if ((norm(
      c(residual_variance, random_effect_variance) -
        c(prev_residual_variance, prev_random_effect_variance),
      "2"
    ) < tol) & (norm(coefs - prev_coefs, "2") < tol)) {
      parameters_converged <- TRUE
    }

    iter <- iter + 1
    iteration_log <- rbind(
      iteration_log,
      c(residual_variance, random_effect_variance, coefs)
    )
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_lmm",
    ## Input/Processed parameters
    cat_init = cat_init,
    tau = tau,
    residual_variance_0 = residual_variance_0,
    random_effect_variance_0 = random_effect_variance_0,
    coefs_0 = coefs_0,
    optimize_domain = optimize_domain,
    max_iter = max_iter,
    tol = tol,
    ## Result
    coefficients = as.numeric(coefs),
    iteration_log = iteration_log,
  )

  class(cat_model) <- c(cat_model$class, "cat_lmm")

  return(cat_model)
}

#' Title
#'
#' @param data_
#' @param tau
#' @param ...
#' @param max_iter
#' @param tol
#' @param init_sig_square
#' @param init_coefficients
#'
#' @return
#' @export
#'
#' @examples
cat_lmm <- function(data_,
                    tau = NULL,
                    ...,
                    max_iter = 500,
                    tol = 1e-08,
                    init_sig_square = c(1, 1),
                    init_coefficients = NULL,
                    optimize_domain = c(0, 30)) {
  .not_true_stop(
    is.data_lmm(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization_lmm."
    )
  )

  tau <- if (is.null(tau)) ncol(data_$obs_x) / 4 else tau

  ZZT <- tcrossprod(as.matrix(data_$obs_z))
  eigen_ZZT <- eigen(ZZT)
  d_seq <- eigen_ZZT$values
  U <- eigen_ZZT$vectors

  Z_synZ_synT <- tcrossprod(as.matrix(data_$syn_z))
  eigen_Z_synZ_synT <- eigen(Z_synZ_synT)
  d_seq_star <- eigen_Z_synZ_synT$values
  U_star <- eigen_Z_synZ_synT$vectors


  iter <- 0
  coefficients_new <- if (is.null(init_coefficients)) rnorm(ncol(data_$obs_x)) else init_coefficients
  sig_square_par_new <- init_sig_square

  parameter_diff <- TRUE

  iter_collect <- matrix(
    data = c(sig_square_par_new, coefficients_new),
    ncol = (2 + length(coefficients_new))
  )
  colnames(iter_collect) <- c("sigma_square", "sigma_u_square", paste0("beta[", 1:length(coefficients_new), "]"))

  while ((iter < max_iter) & (parameter_diff)) {
    coefficients_old <- coefficients_new
    V_inv <- quadform::quad.tform(
      diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq)), U
    )
    V_star_inv <- quadform::quad.tform(
      diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq_star)), U_star
    )
    coefficients_new <- MASS::ginv(tau / data_$syn_size * quadform::quad.form(V_star_inv, data_$syn_x) +
      quadform::quad.form(V_inv, data_$obs_x)) %*% (
      quadform::quad.3form(V_inv, data_$obs_x, data_$obs_y) +
        tau / data_$syn_size * quadform::quad.3form(V_star_inv, data_$syn_x, data_$syn_y))
    q <- t(U) %*% (data_$obs_y - data_$obs_x %*% coefficients_new)
    q_star <- t(U_star) %*% (data_$syn_y - data_$syn_x %*% coefficients_new)
    sig_square_par_old <- sig_square_par_new
    sigma_square_old <- sig_square_par_new[1]
    sigma_u_square_new <- stats::optimize(update_sigma_u_square_cat, optimize_domain, sigma_square_old, d_seq, d_seq_star, q, q_star, tau)$minimum
    sigma_square_new <- stats::optimize(update_sigma_square_cat, optimize_domain, sigma_u_square_new, d_seq, d_seq_star, q, q_star, tau)$minimum
    sig_square_par_new <- c(sigma_square_new, sigma_u_square_new)

    if ((norm(sig_square_par_new - sig_square_par_old, "2") < tol) & (norm(coefficients_new - coefficients_old, "2") < tol)) {
      parameter_diff <- FALSE
    }
    iter <- iter + 1
    iter_collect <- rbind(iter_collect, c(sig_square_par_new, coefficients_new))
  }

  cat_lmer_obj <- list(
    function_name = "cat_lmm",
    data_ = data_,
    var_par = sig_square_par_new,
    coefficients = as.numeric(coefficients_new),
    tau = tau,
    iter_collect = iter_collect,
    max_iter = max_iter,
    tol = tol,
    optimize_domain = optimize_domain
  )

  class(cat_lmer_obj) <- c(cat_lmer_obj$class, "lmm", "cat_lmm")


  return(cat_lmer_obj)
}

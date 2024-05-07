#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param tau
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
surv_glm <-
  function(X,
           Y,
           Z,
           method = c(
             "mle",
             "ridge"
           ),
           ...,
           max_iter = 500,
           tol = 1e-08,
           lambda = 1,
           init_sig_square = c(1, 1),
           init_beta = NULL) {
    # -------------------------------------------------------------------------


    validate_surv_formula(formula)

    p <- ncol(X)
    ZZT <- tcrossprod(Z)
    eigen_ZZT <- eigen(ZZT)
    d_seq <- eigen_ZZT$values
    U <- eigen_ZZT$vectors

    iter <- 0
    beta_new <- if (is.null(init_beta)) rnorm(p) else init_beta
    sig_square_par_new <- init_sig_square
    parameter_diff <- TRUE

    while ((iter < MAX_iter) & (parameter_diff)) {
      beta_old <- beta_new
      V_inv <- quad.tform(diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq)), U)
      beta_new <- switch(method,
        "mle" = ginv(quad.form(V_inv, X)) %*% quad.3form(V_inv, X, Y),
        "ridge" = ginv(lambda * diag(p) + quad.form(V_inv, X)) %*% quad.3form(V_inv, X, Y)
      )
      q <- t(U) %*% (Y - X %*% beta_new)
      sig_square_par_old <- sig_square_par_new
      sigma_square_old <- sig_square_par_new[1]
      sigma_u_square_new <- optimize(update_sigma_u_square, c(0, Inf), sigma_square_old, d_seq, q)$minimum
      sigma_square_new <- optimize(update_sigma_square, c(0, Inf), sigma_u_square_new, d_seq, q)$minimum
      sig_square_par_new <- c(sigma_square_new, sigma_u_square_new)

      if ((norm(sig_square_par_new - sig_square_par_old, "2") < tol) & (norm(beta_new - beta_old, "2") < tol)) {
        parameter_diff <- FALSE
      }
      iter <- iter + 1
    }
    return(list(
      coefficients = as.numeric(beta_new),
      var_par = sig_square_par_new
    ))
  }

# collection of function LMM_ridge, LMM_cat, LMM_MLE


LMM_MLE_Oct4 <- function(X, Y, Z, MAX_iter = 500, tol = 1e-08) {
  p <- ncol(X)
  ZZT <- tcrossprod(Z)
  eigen_ZZT <- eigen(ZZT)
  d_seq <- eigen_ZZT$values
  U <- eigen_ZZT$vectors
  update_sigma_square <- function(sigma_square, sigma_u_square, d_seq, q) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
  }
  update_sigma_u_square <- function(sigma_u_square, sigma_square, d_seq, q) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
  }
  iter <- 0
  beta_new <- rnorm(p)
  sig_square_par_new <- c(1, 1)
  parameter_diff <- TRUE
  while ((iter < MAX_iter) & (parameter_diff)) {
    beta_old <- beta_new
    V_inv <- quad.tform(diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq)), U)
    beta_new <- ginv(quad.form(V_inv, X)) %*% quad.3form(V_inv, X, Y)
    q <- t(U) %*% (Y - X %*% beta_new)
    sig_square_par_old <- sig_square_par_new
    sigma_square_old <- sig_square_par_new[1]
    sigma_u_square_new <- optimize(update_sigma_u_square, c(0, 30), sigma_square_old, d_seq, q)$minimum
    sigma_square_new <- optimize(update_sigma_square, c(0, 30), sigma_u_square_new, d_seq, q)$minimum
    sig_square_par_new <- c(sigma_square_new, sigma_u_square_new)

    if ((norm(sig_square_par_new - sig_square_par_old, "2") < tol) & (norm(beta_new - beta_old, "2") < tol)) {
      parameter_diff <- FALSE
    }
    iter <- iter + 1
  }
  return(list(beta = as.numeric(beta_new), var_par = sig_square_par_new))
}
LMM_ridge_Oct4 <- function(X, Y, Z, lambda, MAX_iter = 500, tol = 1e-08) {
  p <- ncol(X)
  ZZT <- tcrossprod(Z)
  eigen_ZZT <- eigen(ZZT)
  d_seq <- eigen_ZZT$values
  U <- eigen_ZZT$vectors
  update_sigma_square <- function(sigma_square, sigma_u_square, d_seq, q) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
  }
  update_sigma_u_square <- function(sigma_u_square, sigma_square, d_seq, q) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square))
  }
  iter <- 0
  beta_new <- rnorm(p)
  sig_square_par_new <- c(1, 1)
  parameter_diff <- TRUE
  while ((iter < MAX_iter) & (parameter_diff)) {
    beta_old <- beta_new
    V_inv <- quad.tform(diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq)), U)
    beta_new <- ginv(lambda * diag(p) + quad.form(V_inv, X)) %*% quad.3form(V_inv, X, Y)
    q <- t(U) %*% (Y - X %*% beta_new)
    sig_square_par_old <- sig_square_par_new
    sigma_square_old <- sig_square_par_new[1]
    sigma_u_square_new <- optimize(update_sigma_u_square, c(0, 30), sigma_square_old, d_seq, q)$minimum
    sigma_square_new <- optimize(update_sigma_square, c(0, 30), sigma_u_square_new, d_seq, q)$minimum
    sig_square_par_new <- c(sigma_square_new, sigma_u_square_new)

    if ((norm(sig_square_par_new - sig_square_par_old, "2") < tol) & (norm(beta_new - beta_old, "2") < tol)) {
      parameter_diff <- FALSE
    }
    iter <- iter + 1
  }
  return(list(beta = as.numeric(beta_new), var_par = sig_square_par_new))
}
LMM_cat_Oct4 <- function(X, Y, Z, X_syn, Y_syn, Z_syn, tau, MAX_iter = 500, tol = 1e-08) {
  p <- ncol(X)
  M <- nrow(X_syn)

  ZZT <- tcrossprod(Z)
  eigen_ZZT <- eigen(ZZT)
  d_seq <- eigen_ZZT$values
  U <- eigen_ZZT$vectors

  Z_synZ_synT <- tcrossprod(Z_syn)
  eigen_Z_synZ_synT <- eigen(Z_synZ_synT)
  d_seq_star <- eigen_Z_synZ_synT$values
  U_star <- eigen_Z_synZ_synT$vectors



  update_sigma_square <- function(sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) + tau / length(d_seq_star) *
      sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))
  }
  update_sigma_u_square <- function(sigma_u_square, sigma_square, d_seq, d_seq_star, q, q_star, tau) { # minimization program
    sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) + tau / length(d_seq_star) *
      sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))
  }
  iter <- 0
  beta_new <- rnorm(p)
  sig_square_par_new <- c(1, 1)
  parameter_diff <- TRUE
  while ((iter < MAX_iter) & (parameter_diff)) {
    beta_old <- beta_new
    V_inv <- quadform::quad.tform(diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq)), U)
    V_star_inv <- quadform::quad.tform(diag(1 / (sig_square_par_new[1] + sig_square_par_new[2] * d_seq_star)), U_star)
    beta_new <- MASS::ginv(tau / M * quadform::quad.form(V_star_inv, X_syn) +
      quadform::quad.form(V_inv, X)) %*% (
      quadform::quad.3form(V_inv, X, Y) + tau / M * quadform::quad.3form(V_star_inv, X_syn, Y_syn))
    q <- t(U) %*% (Y - X %*% beta_new)
    q_star <- t(U_star) %*% (Y_syn - X_syn %*% beta_new)
    sig_square_par_old <- sig_square_par_new
    sigma_square_old <- sig_square_par_new[1]
    sigma_u_square_new <- stats::optimize(update_sigma_u_square, c(0, 30), sigma_square_old, d_seq, d_seq_star, q, q_star, tau)$minimum
    sigma_square_new <- stats::optimize(update_sigma_square, c(0, 30), sigma_u_square_new, d_seq, d_seq_star, q, q_star, tau)$minimum
    sig_square_par_new <- c(sigma_square_new, sigma_u_square_new)

    if ((norm(sig_square_par_new - sig_square_par_old, "2") < tol) & (norm(beta_new - beta_old, "2") < tol)) {
      parameter_diff <- FALSE
    }
    iter <- iter + 1
  }
  browser()
  return(list(beta = as.numeric(beta_new), var_par = sig_square_par_new))
}

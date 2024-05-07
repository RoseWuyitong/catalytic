# in this script we ensemble all bayesian sampler in one file, only comsider one random effect variance sigma_u_square, diff group same variance


############## version 3 LMM

update_beta_version3_same_u <- function(sigma_square, sigma_u_square, random_effect_list, mydata, syndata, tau) {
  m <- length(mydata$Y_list)
  # n_i_seq<-  sapply(1:m,function(aa){ nrow(mydata$X_list[[aa]])})
  p <- ncol(mydata$X_list[[1]])
  # q<- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])
  V <- matrix(0, p, p)
  b <- matrix(0, nc = 1, nr = p)
  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list

  X_list_syn <- syndata$X_list
  Z_list_syn <- syndata$Z_list
  Y_list_syn <- syndata$Y_list

  for (i in 1:m) {
    V <- V + crossprod(X_list_[[i]]) + tau / M / m * crossprod(X_list_syn[[i]])
    b <- b + t(X_list_[[i]]) %*% (Y_list_[[i]] - Z_list_[[i]] %*% random_effect_list[[i]]) +
      tau / m / M * t(X_list_syn[[i]]) %*% (Y_list_syn[[i]] - Z_list_syn[[i]] %*% random_effect_list[[i]])
  }
  inv_V <- ginv(V)
  betanew <- rmvnorm(1, mean = inv_V %*% b, sigma = sigma_square * inv_V)
  return(as.numeric(betanew))
}
update_random_effect_version3_same_u <- function(beta, q, sigma_u_square, mydata, syndata, ZTZ_list, ZTZ_syn_list, sigma_square, tau, M) {
  m <- length(ZTZ_list)
  lapply(1:m, function(index) {
    Vi <- (1 + tau / m / M) * sigma_square / sigma_u_square * diag(q) + ZTZ_list[[index]] + tau / m / M * ZTZ_syn_list[[index]]
    bi <- t(mydata$Z_list[[index]]) %*% (mydata$Y_list[[index]] - mydata$X_list[[index]] %*% beta) +
      tau / m / M * t(syndata$Z_list[[index]]) %*% (syndata$Y_list[[index]] - syndata$X_list[[index]] %*% beta)
    inv_Vi <- ginv(Vi)
    as.numeric(rmvnorm(1, mean = inv_Vi %*% bi, sigma = inv_Vi * sigma_square))
  })
}

update_sigma_square_version3_same_u <- function(beta, random_effect_list, mydata, syndata, tau, prior_rate = 0) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  # p<-ncol(mydata$X_list[[1]])
  # q<- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])
  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list

  X_list_syn <- syndata$X_list
  Z_list_syn <- syndata$Z_list
  Y_list_syn <- syndata$Y_list
  shape_inv <- 0.5 * (tau + sum(n_i_seq))
  rate_inv <- 0
  for (i in 1:m) {
    rate_inv <- rate_inv + norm(Y_list_[[i]] - X_list_[[i]] %*% beta - Z_list_[[i]] %*% random_effect_list[[i]], "2")^2 +
      tau / m / M * norm(Y_list_syn[[i]] - X_list_syn[[i]] %*% beta - Z_list_syn[[i]] %*% random_effect_list[[i]], "2")^2
  }
  # cat("rate is ",rate_inv/2,"\n")
  # print(c(shape_inv,rate_inv/2))
  rinvgamma(1, shape = shape_inv, rate = prior_rate + rate_inv / 2)
}
update_sigma_u_square_version3_same_u <- function(q, random_effect_list, tau, M, prior_rate = 0) {
  # q=length(random_effect_list[[1]])
  m <- length(random_effect_list)

  rate_inv <- sum(sapply(1:m, function(index) {
    (1 / 2 + 1 / 2 * tau / m / M) * crossprod(random_effect_list[[index]])
  }))
  shapeinv <- m * q / 2 + q / 2 * tau / M - 1
  proposed <- rinvgamma(1, shapeinv, rate = prior_rate + rate_inv)
  # print(c(shapeinv,rate_inv))
  proposed
}

############## version 3 LMM complete




############# traditional prior LMM



update_beta_traditional_prior <- function(sigma_square, sigma_u_square, random_effect_list, mydata) {
  m <- length(mydata$Y_list)
  # n_i_seq<-  sapply(1:m,function(aa){ nrow(mydata$X_list[[aa]])})
  p <- ncol(mydata$X_list[[1]])


  V <- matrix(0, p, p)
  b <- matrix(0, nc = 1, nr = p)
  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list


  for (i in 1:m) {
    V <- V + crossprod(X_list_[[i]])
    b <- b + t(X_list_[[i]]) %*% (Y_list_[[i]] - Z_list_[[i]] %*% random_effect_list[[i]])
  }
  inv_V <- ginv(V)
  betanew <- rmvnorm(1, mean = inv_V %*% b, sigma = sigma_square * inv_V)
  return(as.numeric(betanew))
}
update_random_effect_traditional_prior <- function(beta, q, sigma_u_square, mydata, ZTZ_list, sigma_square) {
  m <- length(ZTZ_list)
  lapply(1:m, function(index) {
    Vi <- sigma_square / sigma_u_square * diag(q) + ZTZ_list[[index]]
    bi <- t(mydata$Z_list[[index]]) %*% (mydata$Y_list[[index]] - mydata$X_list[[index]] %*% beta)
    inv_Vi <- ginv(Vi)
    as.numeric(rmvnorm(1, mean = inv_Vi %*% bi, sigma = inv_Vi * sigma_square))
  })
}

update_sigma_square_traditional_prior <- function(beta, random_effect_list, mydata, prior_rate = 0) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  # p<-ncol(mydata$X_list[[1]])
  # q<- ncol(mydata$Z_list[[1]])

  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list

  shape_inv <- 0.5 * (sum(n_i_seq))
  rate_inv <- 0
  for (i in 1:m) {
    rate_inv <- rate_inv + norm(Y_list_[[i]] - X_list_[[i]] %*% beta - Z_list_[[i]] %*% random_effect_list[[i]], "2")^2
  }

  rinvgamma(1, shape = shape_inv, rate = prior_rate + rate_inv / 2)
}

update_sigma_u_square_traditional_prior <- function(q, random_effect_list, prior_rate = 0) {
  m <- length(random_effect_list)

  rate_inv <- sum(sapply(1:m, function(index) {
    (1 / 2) * crossprod(random_effect_list[[index]])
  }))
  shapeinv <- m * q / 2 - 1
  proposed <- rinvgamma(1, shapeinv, rate = prior_rate + rate_inv)
  # print(c(shapeinv,rate_inv))
  proposed
}


#############

#################### MCMC power prior

update_beta_power_prior <- function(sigma_square, sigma_u_square, random_effect_list,
                                    random_effect_syn_list, mydata, syndata, tau) {
  m <- length(mydata$Y_list)
  # n_i_seq<-  sapply(1:m,function(aa){ nrow(mydata$X_list[[aa]])})
  p <- ncol(mydata$X_list[[1]])
  # q<- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])
  V <- matrix(0, p, p)
  b <- matrix(0, nc = 1, nr = p)
  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list

  X_list_syn <- syndata$X_list
  Z_list_syn <- syndata$Z_list
  Y_list_syn <- syndata$Y_list

  for (i in 1:m) {
    V <- V + crossprod(X_list_[[i]]) + tau / M / m * crossprod(X_list_syn[[i]])
    b <- b + t(X_list_[[i]]) %*% (Y_list_[[i]] - Z_list_[[i]] %*% random_effect_list[[i]]) +
      tau / m / M * t(X_list_syn[[i]]) %*% (Y_list_syn[[i]] - Z_list_syn[[i]] %*% random_effect_syn_list[[i]])
  }
  inv_V <- ginv(V)
  betanew <- rmvnorm(1, mean = inv_V %*% b, sigma = sigma_square * inv_V)
  return(as.numeric(betanew))
}
update_random_effect_power_prior <- function(beta, q, sigma_u_square, mydata, ZTZ_list, sigma_square) {
  m <- length(ZTZ_list)
  lapply(1:m, function(index) {
    Vi <- sigma_square / sigma_u_square * diag(q) + ZTZ_list[[index]]
    bi <- t(mydata$Z_list[[index]]) %*% (mydata$Y_list[[index]] - mydata$X_list[[index]] %*% beta)
    inv_Vi <- ginv(Vi)
    as.numeric(rmvnorm(1, mean = inv_Vi %*% bi, sigma = inv_Vi * sigma_square))
  })
}
update_random_effect_syn_power_prior <- function(beta, q, sigma_u_square, syndata, ZTZ_syn_list, sigma_square, tau, M) {
  m <- length(ZTZ_syn_list)
  lapply(1:m, function(index) {
    Vi <- sigma_square / sigma_u_square * diag(q) + tau / m / M * ZTZ_syn_list[[index]]
    bi <- tau / m / M * t(syndata$Z_list[[index]]) %*% (syndata$Y_list[[index]] - syndata$X_list[[index]] %*% beta)
    inv_Vi <- ginv(Vi)
    as.numeric(rmvnorm(1, mean = inv_Vi %*% bi, sigma = inv_Vi * sigma_square))
  })
}
update_sigma_square_power_prior <- function(beta, random_effect_list,
                                            random_effect_syn_list, mydata, syndata, tau, prior_rate = 0) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  # p<-ncol(mydata$X_list[[1]])
  # q<- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])
  X_list_ <- mydata$X_list
  Z_list_ <- mydata$Z_list
  Y_list_ <- mydata$Y_list

  X_list_syn <- syndata$X_list
  Z_list_syn <- syndata$Z_list
  Y_list_syn <- syndata$Y_list
  shape_inv <- 0.5 * (tau + sum(n_i_seq))
  rate_inv <- 0
  for (i in 1:m) {
    rate_inv <- rate_inv + norm(Y_list_[[i]] - X_list_[[i]] %*% beta - Z_list_[[i]] %*% random_effect_list[[i]], "2")^2 +
      tau / m / M * norm(Y_list_syn[[i]] - X_list_syn[[i]] %*% beta - Z_list_syn[[i]] %*% random_effect_syn_list[[i]], "2")^2
  }
  # cat("rate is ",rate_inv/2,"\n")
  # print(c(shape_inv,rate_inv/2))
  rinvgamma(1, shape = shape_inv, rate = prior_rate + rate_inv / 2)
}
update_sigma_u_square_power_prior <- function(q, random_effect_list,
                                              random_effect_syn_list, tau, M, prior_rate = 0) {
  # q=length(random_effect_list[[1]])
  m <- length(random_effect_list)

  rate_inv <- sum(sapply(1:m, function(index) {
    (crossprod(random_effect_syn_list[[index]]) / 2 + crossprod(random_effect_list[[index]]) / 2)
  }))
  shapeinv <- m * q / 2 + m * q / 2 - 1
  proposed <- rinvgamma(1, shapeinv, rate = rate_inv + prior_rate)
  # print(c(shapeinv,rate_inv))
  proposed
}
####################



########  no random effect, aug5, version 2

log_posterior_ver2_no_u <- function(beta, sigma_square, sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  p <- ncol(mydata$X_list[[1]])
  q <- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])

  A_list <- lapply(1:m, function(aa) {
    sigma_u_square * ZZT_list[[aa]] + sigma_square * diag(n_i_seq[aa])
  })

  A_list_star <- lapply(1:m, function(aa) {
    sigma_u_square * ZZT_syn_list[[aa]] + sigma_square * diag(M)
  })

  sum(sapply(1:m, function(i) {
    -1 / 2 * determinant(A_list[[i]])$modulus - 1 / 2 * tau / m / M * determinant(A_list_star[[i]])$modulus - 1 / 2 * quad.form.inv(A_list[[i]], mydata$Y_list[[i]] - mydata$X_list[[i]] %*% beta) -
      1 / 2 * quad.form.inv(A_list_star[[i]], syndata$Y_list[[i]] - syndata$X_list[[i]] %*% beta) * tau / m / M - log(sigma_u_square)
  })) - log(sigma_square)
}


update_beta_ver2_no_u <- function(beta, sigma_square, sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  p <- ncol(mydata$X_list[[1]])
  q <- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])

  beta_old <- beta

  A_inv_list <- lapply(1:m, function(aa) {
    ginv(sigma_u_square * ZZT_list[[aa]] + sigma_square * diag(n_i_seq[aa]))
  })

  A_inv_list_star <- lapply(1:m, function(aa) {
    ginv(sigma_u_square * ZZT_syn_list[[aa]] + sigma_square * diag(M))
  })
  V <- matrix(0, p, p)
  V_syn <- V

  b <- matrix(0, nc = 1, nr = p)
  b_syn <- b

  for (index in 1:m) {
    V <- V + quad.form(A_inv_list[[index]], mydata$X_list[[index]])
    V_syn <- V_syn + quad.form(A_inv_list_star[[index]], syndata$X_list[[index]])
    b <- b + quad.3form(A_inv_list[[index]], mydata$X_list[[index]], mydata$Y_list[[index]])
    b_syn <- b_syn + quad.3form(A_inv_list_star[[index]], syndata$X_list[[index]], syndata$Y_list[[index]])
  }

  cov_beta <- ginv(V + tau / m / M * V_syn)
  b_combine <- b + tau / m / M * b_syn
  beta <- as.numeric(rmvnorm(1, mean = cov_beta %*% b_combine, sigma = cov_beta))

  return(list(beta_old = beta_old, beta_new = beta))
}

update_sigma_square_ver2_no_u <- function(beta, sigma_square, sigma_u_square, mydata, tau, syndata, MH_sd = 0.1, ZZT_list, ZZT_syn_list) {
  sigma_square_old <- sigma_square
  cur_value <- log(sigma_square_old)
  # we will use random walk MH to sample once
  prop_value <- rnorm(1, mean = cur_value, sd = MH_sd)
  prop_sigma_square <- exp(prop_value)

  log_posterior_cur <- log_posterior_ver2_no_u(beta, sigma_square_old, sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list)
  log_posterior_prop <- log_posterior_ver2_no_u(beta, prop_sigma_square, sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list)

  logR <- log_posterior_prop - log_posterior_cur

  u <- log(runif(1)) < logR

  sigma_square_new <- sigma_square_old
  if (u == 1) {
    sigma_square_new <- prop_sigma_square
  }

  return(list(sigma_square_new = sigma_square_new, accept = u))
}


update_sigma_u_square_ver2_no_u <- function(beta, sigma_square, sigma_u_square, mydata, tau, syndata, MH_sd = 0.1, ZZT_list, ZZT_syn_list) {
  m <- length(mydata$Y_list)
  n_i_seq <- sapply(1:m, function(aa) {
    nrow(mydata$X_list[[aa]])
  })
  p <- ncol(mydata$X_list[[1]])
  q <- ncol(mydata$Z_list[[1]])
  M <- nrow(syndata$X_list[[1]])
  accept_ <- 0

  cur_sigma_u_square <- sigma_u_square
  cur_value <- log(cur_sigma_u_square)

  prop_value <- rnorm(1, mean = cur_value, sd = MH_sd)
  prop_sigma_u_square <- exp(prop_value)

  log_posterior_cur <- log_posterior_ver2_no_u(beta, sigma_square, cur_sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list)
  log_posterior_prop <- log_posterior_ver2_no_u(beta, sigma_square, prop_sigma_u_square, mydata, tau, syndata, ZZT_list, ZZT_syn_list)


  logR <- log_posterior_prop - log_posterior_cur

  u <- log(runif(1)) < logR

  alpha_i_sigma_square_new <- cur_sigma_u_square
  if (u == 1) {
    sigma_u_square <- prop_sigma_u_square
  }
  accept_ <- u


  return(list(sigma_u_square = sigma_u_square, accept_seq = accept_))
}
########

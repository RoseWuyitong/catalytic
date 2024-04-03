#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param ...
#' @param iter
#' @param init_tau
#' @param init_beta
#' @param gibbs_iter
#' @param algorithm
#' @param alpha
#' @param gamma
#' @param stan_refresh
#' @param warmup
#'
#' @return
#' @export
#'
#' @examples
cat_gibbs_glm <- function(formula,
                          family = gaussian(),
                          data_,
                          init_tau = NULL,
                          init_beta = NULL,
                          ...,
                          iter = 5,
                          warmup = 4,
                          gibbs_iter = 100,
                          algorithm = "HMC",
                          alpha = 2,
                          gamma = 1,
                          stan_refresh = 0) {
  .not_in_stop(
    c(
      "obs_size", "obs_data", "obs_x", "obs_y", "syn_size", "syn_data", "syn_x",
      "syn_y", "syn_model", "data", "size", "x", "y"
    ),
    names(data_),
    "Not all necessary variables are present in data_."
  )

  data_$adj_x <- model.frame(rhs(formula, TRUE), data_$x)
  data_$adj_obs_x <- data_$adj_x[1:data_$obs_size, ]
  data_$adj_syn_x <- data_$adj_x[(data_$obs_size + 1):data_$size, ]

  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- get_family(family)

  kappa_ <- log_density(fam)(
    data_$adj_syn_x,
    data_$syn_y,
    stats::glm(formula = formula, data = data_$syn_data, family = family),
    weights = 1 / data_$syn_size)

  data <- list(
    obs_n = data_$obs_size,
    syn_n = data_$syn_size,
    x_col_n = ncol(data_$adj_x),
    tau = if (is.null(init_tau)) ncol(data_$adj_x) / 4 else init_tau,
    obs_x = data_$adj_obs_x,
    obs_y = as.numeric(data_$obs_y[[1]]),
    syn_x = data_$adj_syn_x,
    syn_y = as.numeric(data_$syn_y[[1]]),
    kappa = kappa_
  )

  stan_model <- stan_model_selection(
    fam = fam,
    joint = TRUE,
    gibbs = TRUE,
    use_theta = FALSE
  )

  if (is.null(init_beta)) {
    sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      chains = 1,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm,
      refresh = stan_refresh
    )

    init_beta <- rstan::extract(sample_model)$beta[
      nrow(rstan::extract(sample_model)$beta),
    ]
  }

  beta_collect <- c(init_beta)
  tau_collect <- c(data$tau)

  for (gibb in 1:gibbs_iter) {
    cat("\ngibbs: ", gibb)

    sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      init = function(...) list(beta = init_beta),
      chains = 1,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm,
      refresh=stan_refresh
    )

    init_beta <- rstan::extract(sample_model)$beta[
      nrow(rstan::extract(sample_model)$beta),
    ]

    beta_collect <- append(beta_collect, init_beta)

    data$tau <- rgamma(
      n = 1,
      shape = ncol(data_$adj_x) + alpha,
      scale = kappa_ +
        1 / gamma -
        log_density(fam)(
          data_$adj_syn_x,
          data_$syn_y,
          beta = init_beta,
          weights = 1 / data_$syn_size)
    )

    tau_collect <- append(tau_collect, data$tau)
  }

  cat_gibbs <- list(
    cat_function = "cat_gibbs_glm",
    formula = formula,
    data_ = data_,
    tau = data$tau,
    family = family,
    joint = TRUE,
    algorithm = algorithm,
    sample_model = sample_model,
    beta_collect = matrix(beta_collect,
      ncol = (ncol(data_$adj_x) + 1),
      byrow = TRUE
    ),
    tau_collect = tau_collect,
    beta = init_beta,
    iter = iter,
    warmup = warmup,
    alpha=alpha,
    gamma=gamma,
    gibbs_iter = gibbs_iter
  )

  class(cat_gibbs) <- c(cat_gibbs$class, "cat_gibbs_glm", "cat_glm")

  return(cat_gibbs)
}


#' Title
#'
#' @param cat_gibbs
#' @param gibbs_iter
#'
#' @return
#' @export
#'
#' @examples
add_gibbs_iter.cat_gibbs_glm <- function(cat_gibbs, gibbs_iter=100) {
  updated_cat_gibbs <- cat_gibbs_glm(formula = cat_gibbs$formula,
                                      family = cat_gibbs$family,
                                      data_ = cat_gibbs$data_,
                                      init_tau = cat_gibbs$tau,
                                      init_beta = cat_gibbs$beta,
                                      iter = cat_gibbs$iter,
                                      warmup = cat_gibbs$warmup,
                                      gibbs_iter = gibbs_iter,
                                      algorithm = cat_gibbs$algorithm,
                                      alpha = cat_gibbs$alpha,
                                      gamma = cat_gibbs$gamma)

  updated_cat_gibbs$gibbs_iter <-  cat_gibbs$gibbs_iter + gibbs_iter
  updated_cat_gibbs$beta_collect <- rbind(cat_gibbs$beta_collect, updated_cat_gibbs$beta_collect[-1,])
  updated_cat_gibbs$tau_collect <- c(cat_gibbs$tau_collect, updated_cat_gibbs$tau_collect[-1])
  return(updated_cat_gibbs)
}

add_gibbs_iter <- function(cat_gibbs, gibbs_iter) {
  UseMethod("add_gibbs_iter")
}

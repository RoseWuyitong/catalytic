#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param tau
#' @param joint
#' @param ...
#' @param chains
#' @param iter
#' @param warmup
#' @param mcmc_algorithm
#' @param joint_alpha
#' @param joint_gamma
#'
#' @return
#' @export
#'
#' @examples
cat_bayesian_glm <- function(formula,
                             family = gaussian(),
                             data_,
                             tau = NULL,
                             joint = FALSE,
                             gibbs = FALSE,
                             use_theta = FALSE,
                             ...,
                             chains = 4,
                             iter = 2000,
                             joint_lower = 0.05,
                             warmup = 1000,
                             algorithm = "HMC",
                             joint_alpha = 2,
                             joint_gamma = 1,
                             gibbs_iter = 1000) {
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


  if (is.null(tau)) tau <- ncol(data_$adj_x) / 4

  kappa_ <- log_density(fam)(
    data_$adj_syn_x,
    data_$syn_y,
    stats::glm(formula = formula, data = data_$syn_data, family = family),
    weights = 1 / data_$syn_size)

  data <- list(
    obs_n = data_$obs_size,
    syn_n = data_$syn_size,
    x_col_n = ncol(data_$adj_x),
    tau = if (is.null(tau)) ncol(data_$adj_x) / 4 else tau,
    obs_x = data_$adj_obs_x,
    obs_y = as.numeric(data_$obs_y[[1]]),
    syn_x = data_$adj_syn_x,
    syn_y = as.numeric(data_$syn_y[[1]]),
    lower_ = joint_lower,
    kappa = kappa_,
    alpha = joint_alpha,
    gamma = joint_gamma
  )


  stan_model <- stan_model_selection(fam = fam,
                                     joint = joint,
                                     gibbs = gibbs,
                                     use_theta = use_theta)

  if (fam == "binomial" & joint & gibbs) {
    sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm,
    )

    last_beta <- rstan::extract(sample_model)$beta[
      nrow(rstan::extract(sample_model)$beta),
    ]

    beta_collect <- c(last_beta)
    tau_collect <- c(tau)

    for (gibb in 1:gibbs_iter) {
      cat("\ngibbs: ", gibb)

      sample_model <- rstan::sampling(
        object = stan_model,
        data = data,
        init = function(...) list(beta = last_beta),
        chains = chains,
        iter = iter,
        warmup = warmup,
        algorithm = algorithm
      )

      last_beta <- rstan::extract(sample_model)$beta[
        nrow(rstan::extract(sample_model)$beta),
      ]

      beta_collect <- append(beta_collect, last_beta)

      data$tau <- rgamma(
        n = 1,
        shape = ncol(data_$adj_x) + joint_alpha,
        scale = kappa_ +
          1 / joint_gamma -
          log_density(fam)(
            data_$adj_syn_x,
            data_$syn_y,
            beta = last_beta,
            weights = 1 / data_$syn_size)
      )

      tau_collect <- append(tau_collect, data$tau)
    }
  } else {
    sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm
    )
  }

  cat_bayesian <- list(
    cat_function = "cat_bayesian_glm",
    formula = formula,
    data_ = data_,
    tau = tau,
    family = family,
    joint = joint,
    algorithm = algorithm,
    sample_model = sample_model
  )

  if (fam == "binomial" & joint & gibbs) {
    cat_bayesian <- append(
      list(
        beta_collect = matrix(beta_collect,
          ncol = (ncol(data_$adj_x) + 1),
          byrow = TRUE
        ),
        tau_collect = tau_collect,
        beta = last_beta
      ),
      cat_bayesian
    )
  } else {
    append(list(
      beta = rstan::summary(sample_model)$summary[
        grep("beta", rownames(rstan::summary(sample_model)$summary)),
        "mean"
      ]
    ), cat_bayesian)
  }

  class(cat_bayesian) <- c(cat_bayesian$class, "cat_bayesian_glm", "cat_glm")

  return(cat_bayesian)
}



stan_model_selection <- memoise::memoise(function(fam,
                        joint,
                        gibbs,
                        use_theta) {

  return(rstan::stan_model(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      fam,
      if (joint & isFALSE(gibbs)) "_joint",
      "_bayesian",
      if (fam == "binomial" & joint & use_theta & isFALSE(gibbs)) "_theta",
      if (fam == "binomial" & joint & isFALSE(use_theta) & isFALSE(gibbs)) "_tau",
      ".stan"
    )))
})

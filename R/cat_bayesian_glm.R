cat_bayesian_glm <- function(formula,
                             family = gaussian(),
                             data_,
                             tau = NULL,
                             joint = FALSE,
                             ...,
                             chains = 4,
                             iter = 2000,
                             warmup = 1000,
                             mcmc_algorithm = "HMC",
                             joint_alpha = 2,
                             joint_gamma = 1) {
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
  fam <- family$family

  if (is.null(tau) & !joint) tau <- ncol(data_$x) / 4

  call <- match.call(expand.dots = TRUE)

  rstan_model <- rstan::stan(
    file = paste0(
      "~/catalytic/R/rstan_script/",
      fam,
      if (joint) "_joint",
      "_bayesian.stan"
    ),
    data = list(
      obs_n = data_$obs_size,
      syn_n = data_$syn_size,
      x_col_n = ncol(data_$adj_x),
      obs_x = data_$adj_obs_x,
      obs_y = as.numeric(data_$obs_y[[1]]),
      syn_x = data_$adj_syn_x,
      syn_y = as.numeric(data_$syn_y[[1]]),
      tau = tau,
      kappa = log_density(fam)(data_$adj_syn_x,
                               data_$syn_y,
                               stats::glm(formula=formula, data=data_$syn_data, family = family),
                               weights = 1 / data_$syn_size),
      alpha = joint_alpha,
      gamma = joint_gamma
    ),
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = mcmc_algorithm
  )
  return(rstan_model)
  call["data_"] <- NULL
  call["tau"] <- NULL
  call["family"] <- NULL
  call["joint"] <- NULL
  cat_bayesian =  c(list(data_ = data_,
                         tau = tau,
                         family = family,
                         joint = joint,
                         arg_mean = rstan::summary(rstan_model)$summary[,"mean"],
                         rstan_model = rstan_model), as.list(call[-1]))

  class(cat_bayesian) <- c(cat_bayesian$class, "cat_bayesian_glm")

  return(cat_bayesian)
}

print.cat_bayesian_glm <- function(x, ...) {
  cat("\nFormula for tau point estimate:  ")
  print(x$formula)
  cat(
    "\nMCMC algorithm: ", x$mcmc_algorithm,
    "\nJoint tau: ", x$joint,
    "\n\nMean of estimated argments: ", x$arg_mean
  )
  invisible(x)
}

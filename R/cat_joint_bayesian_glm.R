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
#' @param use_theta
#' @param algorithm
#' @param alpha
#' @param gamma
#' @param tau_lower
#' @param warmup
#'
#' @return
#' @export
#'
#' @examples
cat_joint_bayesian_glm <- function(formula,
                                   family = gaussian(),
                                   data_,
                                   tau = NULL,
                                   ...,
                                   chains = 4,
                                   iter = 2000,
                                   warmup = 1000,
                                   algorithm = "HMC",
                                   alpha = 2,
                                   gamma = 1,
                                   use_theta = FALSE,
                                   tau_lower = 0.05,
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
    obs_x = data_$adj_obs_x,
    obs_y = as.numeric(data_$obs_y[[1]]),
    syn_x = data_$adj_syn_x,
    syn_y = as.numeric(data_$syn_y[[1]]),
    lower_ = tau_lower,
    kappa = kappa_,
    alpha = alpha,
    gamma = gamma
  )


  stan_model <- stan_model_selection(
    fam = fam,
    joint = TRUE,
    gibbs = FALSE,
    use_theta = use_theta
  )

  sample_model <- rstan::sampling(
    object = stan_model,
    data = data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    refresh = stan_refresh
  )

  # TODO: add analysis and suggestion

  cat_joint_bayesian <- list(
    cat_function = "cat_joint_bayesian_glm",
    formula = formula,
    data_ = data_,
    tau = if (use_theta) {
      1 / rstan::summary(sample_model)$summary[
        "theta", "mean"
      ]
    } else {
      rstan::summary(sample_model)$summary[
        "tau", "mean"
      ]
    },
    family = family,
    algorithm = algorithm,
    alpha = alpha,
    gamma = gamma,
    sample_model = sample_model,
    beta = rstan::summary(sample_model)$summary[
      grep("beta", rownames(rstan::summary(sample_model)$summary)),
      "mean"
    ]
  )
  class(cat_joint_bayesian) <- c(cat_joint_bayesian$class, "cat_joint_bayesian_glm", "cat_glm")

  return(cat_joint_bayesian)
}

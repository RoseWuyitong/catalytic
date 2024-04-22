#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param ...
#' @param chains
#' @param iter
#' @param use_theta
#' @param algorithm
#' @param alpha
#' @param gamma
#' @param tau_lower
#' @param warmup
#' @param stan_refresh
#'
#' @return
#' @export
#'
#' @examples
cat_glm_bayes_joint <- function(formula,
                                family = gaussian(),
                                data_,
                                chains = 4,
                                iter = 2000,
                                warmup = 1000,
                                algorithm = "HMC",
                                stan_model_file = NULL,
                                ...,
                                alpha = 2,
                                gamma = 1,
                                use_theta = FALSE,
                                tau_lower = 0.05,
                                stan_refresh = 0) {
  .not_true_stop(
    is.data_(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization."
    )
  )

  data_$adj_x <- model.frame(rhs(formula, TRUE), data_$x)
  data_$adj_obs_x <- data_$adj_x[1:data_$obs_size, ]
  data_$adj_syn_x <- data_$adj_x[(data_$obs_size + 1):data_$size, ]

  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- get_family(family)

  kappa_ <- get_log_density(fam)(
    X = data_$adj_syn_x,
    Y = data_$syn_y,
    model = stats::glm(formula = formula, data = data_$syn_data, family = family),
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

  stan_model <- if (is.null(stan_model_file)) stan_model_selection(
    fam = fam,
    joint = TRUE,
    gibbs = FALSE,
    use_theta = use_theta
  ) else rstan::stan_model(file = stan_model_file)

  sample_model <- rstan::sampling(
    object = stan_model,
    data = data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    refresh = stan_refresh
  )


  if (fam == "binomial") {
    if (use_theta) {
      theta_inform <- rstan::summary(sample_model)$summary["theta", ]

      est_params <- stats::optim(
        par = c(1, 1),
        fn = function(params, percentiles, mean_) {
          alpha <- exp(params[1])
          beta <- exp(params[2])
          calc_percentiles <- invgamma::qinvgamma(c(0.25, 0.5, 0.75),
            shape = alpha,
            scale = beta
          )
          error <- sum((percentiles - calc_percentiles)^2) +
            (mean_ - (beta / (alpha - 1)))^2
          return(error)
        },
        percentiles = c(
          theta_inform["25%"],
          theta_inform["50%"],
          theta_inform["75%"]
        ),
        mean_ = theta_inform["mean"]
      )

      if (invgamma::qinvgamma(c(0.975),
        shape = exp(est_params$par[1]),
        scale = exp(est_params$par[2])
      ) < theta_inform["97.5%"]) {
        if (any(rstan::summary(sample_model)$summary[, "Rhat"] > 1.1)) {
          warning("SUGGESTION: consider using more iteration (e.g.", iter + 1000, ")")
        } else if (alpha < 3) {
          warning("SUGGESTION: consider using larger alpha (e.g.", min(alpha + 0.1, 3), ")")
        }
      }
    } else {
      tau_inform <- rstan::summary(sample_model)$summary["tau", ]

      est_params <- stats::optim(
        par = c(1, 1),
        fn = function(params, percentiles, mean_) {
          alpha <- exp(params[1])
          beta <- exp(params[2])
          calc_percentiles <- stats::qgamma(c(0.25, 0.5, 0.75),
            shape = exp(params[1]),
            scale = exp(params[2])
          )
          error <- sum((percentiles - calc_percentiles)^2) +
            (mean_ - (alpha / beta))^2
          return(error)
        },
        percentiles = c(
          tau_inform["25%"],
          tau_inform["50%"],
          tau_inform["75%"]
        ),
        mean_ = tau_inform["mean"]
      )

      if (stats::qgamma(c(0.025),
        shape = exp(est_params$par[1]),
        scale = exp(est_params$par[2])
      ) > tau_inform["2.5%"]) {
        if (any(rstan::summary(sample_model)$summary[, "Rhat"] > 1.2)) {
          warning("SUGGESTION: Please consider using more iteration (e.g.", iter + 1000, ")")
        } else if (tau_lower < 1) {
          warning("SUGGESTION: Please consider using larger tau_lower (e.g.", min(tau_lower + 0.05, 1), ")")
        } else {
          warning("SUGGESTION: Please consider using theta method.")
        }
      }
    }
  }

  cat_glm_obj <- list(
    function_name = "cat_glm_bayes_joint",
    formula = formula,
    data_ = data_,
    family = family,
    algorithm = algorithm,
    tau = if (use_theta) {
      1 / rstan::summary(sample_model)$summary[
        "theta", "mean"
      ]
    } else {
      rstan::summary(sample_model)$summary[
        "tau", "mean"
      ]
    },
    alpha = alpha,
    gamma = gamma,
    sample_model = sample_model,
    coefficients = rstan::summary(sample_model)$summary[
      grep("beta", rownames(rstan::summary(sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_glm_obj) <- c(cat_glm_obj$class, "bayes", "cat_glm_bayes_joint")

  return(cat_glm_obj)
}

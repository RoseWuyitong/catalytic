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
                             ...,
                             chains = 4,
                             iter = 2000,
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
    kappa = kappa_,
    alpha = joint_alpha,
    gamma = joint_gamma
  )

  if (fam == "binomial" & joint & gibbs) {
    rstan_model <- rstan::stan(
      file = paste0(
        "~/catalytic/R/rstan_script/binomial_bayesian.stan"
      ),
      data = data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm
    )

    last_beta <- rstan::extract(rstan_model)$beta[
      nrow(rstan::extract(rstan_model)$beta),
    ]

    beta_collect <- c(last_beta)
    tau_collect <- c(tau)

    for (gibb in 1:gibbs_iter) {
      cat("\n\ngibbs: ", gibb)

      rstan_model <- rstan::stan(
        file = paste0(
          "~/catalytic/R/rstan_script/binomial_bayesian.stan"
        ),
        data = data,
        chains = chains,
        iter = iter,
        warmup = warmup,
        algorithm = algorithm,
        init = function(...) {
          list(beta = last_beta)
        }
      )


      last_beta <- rstan::extract(rstan_model)$beta[
        nrow(rstan::extract(rstan_model)$beta),
      ]
      cat("last_beta:", last_beta)
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
    rstan_model <- rstan::stan(
      file = paste0(
        "~/catalytic/R/rstan_script/",
        fam,
        if (joint) "_joint",
        "_bayesian.stan"
      ),
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
    rstan_model = rstan_model
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
      beta = rstan::summary(rstan_model)$summary[
        grep("beta", rownames(rstan::summary(rstan_model)$summary)),
        "mean"
      ]
    ), cat_bayesian)
  }
  class(cat_bayesian) <- c(cat_bayesian$class, "cat_bayesian_glm", "cat_glm")

  return(cat_bayesian)
}

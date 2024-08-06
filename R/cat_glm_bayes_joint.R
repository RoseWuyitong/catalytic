#' Bayesian Generalized Linear Model with Joint Tau
#'
#' Fits a Bayesian Generalized Linear Model (GLM) using a joint tau approach with Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @param formula The model formula specifying the relationship between the response variable and predictors.
#'   Use format `RESPONSE ~ COVARIATES`. Alternatively, `~ COVARIATES` (e.g., `~ .`) can be used without specifying the response name.
#' @param family The GLM family to use for modeling. If not provided (`NULL`), defaults to the family specified in `cat_init`.
#' @param cat_init Initialization object created using `cat_glm_initialization`, containing observed and synthetic data, plus other relevant information.
#' @param ... Additional arguments passed to the `rstan::sampling` function for customizing the MCMC sampling process.
#' @param stan_chains Number of Markov chains to run during MCMC sampling. Default is 4.
#' @param stan_iter Total number of iterations per chain in the MCMC sampling process. Default is 2000.
#' @param stan_warmup Number of warmup iterations in the MCMC sampling process, discarded as burn-in. Default is 1000.
#' @param stan_algorithm Sampling algorithm used in `rstan::sampling`. Default is "NUTS" (No-U-Turn Sampler).
#' @param alpha Alpha parameter controlling degrees of freedom for the Student's t-distribution in the joint tau approach. Default is 2.
#' @param gamma Gamma parameter in the joint tau approach. Default is 1.
#' @param binomial_use_theta Logical. If TRUE, use theta parameter in the binomial model. Default is FALSE.
#' @param binomial_joint_alpha Logical. If TRUE, use joint alpha in the binomial model. Default is FALSE.
#' @param binomial_tau_lower Lower limit for the tau parameter in the binomial model. Default is 0.05.
#'
#' @return A list containing various attributes:
#' \itemize{
#'   \item \code{function_name}: Name of the function (`cat_glm_bayes_joint`).
#'   \item \code{formula}: Full formula used for modeling.
#'   \item \code{family}: GLM family used for modeling.
#'   \item \code{cat_init}: Initialization object created from `cat_glm_initialization`.
#'   \item \code{stan_chains}: Number of MCMC chains.
#'   \item \code{stan_iter}: Total iterations per chain.
#'   \item \code{stan_warmup}: Warmup iterations.
#'   \item \code{stan_algorithm}: Sampling algorithm used.
#'   \item \code{alpha}: Alpha parameter for joint tau approach.
#'   \item \code{gamma}: Gamma parameter for joint tau approach.
#'   \item \code{binomial_use_theta}: Whether theta parameter is used in binomial model.
#'   \item \code{binomial_joint_alpha}: Whether joint alpha is used in binomial model.
#'   \item \code{binomial_tau_lower}: Lower limit for tau in binomial model.
#'   \item \code{tau}: Estimated tau parameter.
#'   \item \code{kappa_}: Calculated kappa value.
#'   \item \code{stan_model}: Stan model object.
#'   \item \code{stan_sample_model}: Result object from MCMC sampling.
#'   \item \code{coefficients}: Estimated coefficients from Bayesian GLM model.
#' }
#'
#' @export
#'
#' @examples
#' @seealso
#' \code{\link{cat_glm_initialization}}
cat_glm_bayes_joint <- function(
    formula,
    family = NULL,
    cat_init,
    ...,
    stan_chains = 4,
    stan_iter = 2000,
    stan_warmup = 1000,
    stan_algorithm = "NUTS",
    alpha = 2,
    gamma = 1,
    binomial_use_theta = FALSE,
    binomial_joint_alpha = FALSE,
    binomial_tau_lower = 0.05) {
  if (!is.cat_glm_init(cat_init)) {
    stop(
      "Please enter `cat_init` that generated from `catalytic::cat_glm_initialization`."
    )
  }

  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  f_lhs <- get_formula_lhs(formula)
  if (f_lhs != "NULL" & f_lhs != colnames(cat_init$y)) stop("Please enter the correct y name for formula.")

  cat_init <- get_adjusted_data(cat_init = cat_init, formula = formula)

  if (is.null(family)) family <- cat_init$family

  full_formula <- as.formula(paste0(colnames(cat_init$y), f_rhs))

  suppressWarnings(
    kappa_ <- get_glm_log_density(
      family = family,
      X = cat_init$adj_syn_x,
      Y = cat_init$syn_y,
      coefs = coef(do.call(
        stats::glm,
        list(
          formula = full_formula,
          family = family,
          data = cat_init$syn_data
        )
      )),
      weights = 1 / cat_init$syn_size
    )
  )

  stan_data <- list(
    obs_n = cat_init$obs_size,
    syn_n = cat_init$syn_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_x = cat_init$adj_syn_x,
    obs_y = as.vector(cat_init$obs_y),
    syn_y = as.vector(cat_init$syn_y),
    tau_lower = binomial_tau_lower,
    kappa = kappa_,
    alpha = alpha,
    gamma = gamma
  )

  stan_model <- get_stan_model(
    type = "glm",
    glm_family = family,
    joint_tau = TRUE,
    glm_binomial_joint_theta = binomial_use_theta,
    glm_binomial_joint_alpha = binomial_joint_alpha
  )

  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    algorithm = stan_algorithm,
    ...
  )

  give_glm_bayes_joint_suggestion(
    family,
    stan_iter,
    alpha,
    stan_sample_model,
    binomial_use_theta,
    binomial_joint_alpha,
    binomial_tau_lower
  )

  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]

  if (!is.null(coefs)) {
    names(coefs) <- c("(Intercept)", colnames(cat_init$adj_x))
  }

  cat_glm_bayes_joint_ <- list(
    function_name = "cat_glm_bayes_joint",
    formula = full_formula,
    family = family,
    cat_init = cat_init,
    stan_chains = stan_chains,
    stan_iter = stan_iter,
    stan_warmup = stan_warmup,
    stan_algorithm = stan_algorithm,
    alpha = alpha,
    gamma = gamma,
    binomial_use_theta = binomial_use_theta,
    binomial_joint_alpha = binomial_joint_alpha,
    binomial_tau_lower = binomial_tau_lower,
    tau = if (binomial_use_theta) {
      1 / rstan::summary(stan_sample_model)$summary[
        "theta", "mean"
      ]
    } else {
      rstan::summary(stan_sample_model)$summary[
        "tau", "mean"
      ]
    },
    kappa_ = kappa_,
    stan_model = stan_model,
    stan_sample_model = stan_sample_model,
    coefficients = coefs,
    ...
  )

  class(cat_glm_bayes_joint_) <- c(cat_glm_bayes_joint_$class, "cat_glm_bayes_joint", "cat_glm_bayes", "bayes")

  return(cat_glm_bayes_joint_)
}

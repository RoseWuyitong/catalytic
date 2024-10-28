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
#' @param binomial_joint_theta Logical. If TRUE, use theta parameter in the binomial model. Default is FALSE.
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
#'   \item \code{binomial_joint_theta}: Whether theta parameter is used in binomial model.
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
    cat_init,
    binomial_tau_lower = 0.05,
    binomial_joint_theta = FALSE,
    binomial_joint_alpha = FALSE,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    algorithm = "NUTS",
    tau_alpha = 2,
    tau_gamma = 1,
    gaussian_variance_alpha = NULL,
    gaussian_variance_beta = NULL,
    ...) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    binomial_tau_lower = binomial_tau_lower,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  stan_data <- c(
    list(
      obs_n = cat_init$obs_size,
      syn_n = cat_init$syn_size,
      col_n = ncol(cat_init$adj_x),
      obs_x = cat_init$adj_obs_x,
      syn_x = cat_init$adj_syn_x,
      obs_y = cat_init$obs_y,
      syn_y = cat_init$syn_y,
      kappa = get_glm_log_density(
        family_string = cat_init$family,
        X = cat_init$adj_syn_x,
        Y = cat_init$syn_y,
        coefs = coef(
          do.call(
            stats::glm,
            list(
              formula = full_formula,
              family = cat_init$family,
              data = cat_init$syn_data
            )
          )
        ),
        weights = 1 / cat_init$syn_size
      ),
      tau_alpha = tau_alpha,
      tau_gamma = tau_gamma
    ),
    if (cat_init$family == "gaussian") {
      # Include Gaussian-specific parameters if the family is Gaussian
      if (is.null(cat_init$custom_variance)) {
        # Include the alpha and beta for the inverse-gamma prior for unknown variance
        list(
          variance_alpha = ifelse(is.null(gaussian_variance_alpha), ncol(cat_init$adj_x), gaussian_variance_alpha),
          variance_beta = ifelse(is.null(gaussian_variance_beta), ncol(cat_init$adj_x), gaussian_variance_beta)
        )
      } else {
        # Include sigma value if user set the variance is known
        list(
          sigma = sqrt(cat_init$custom_variance)
        )
      }
      # Include Gaussian-specific parameters if the family is Gaussian
      # if (cat_init$gaussian_known_variance) {
      #   # Include sigma value if user set the variance is known
      #   list(
      #     sigma = sqrt(
      #       ifelse(
      #         is.null(cat_init$custom_variance),
      #         get_glm_custom_var(
      #           cat_init = cat_init,
      #           tau = ncol(cat_init$adj_x), # Specifying tau as the number of predictors for variance adjustment
      #           formula = full_formula
      #         ),
      #         cat_init$custom_variance
      #       )
      #     )
      #   )
      # } else {
      #   # Include the alpha and beta for the inverse-gamma prior for unknown variance
      #   list(
      #     variance_alpha = ncol(cat_init$adj_x),
      #     variance_beta = ncol(cat_init$adj_x)
      #   )
      # }
    } else {
      # Include Binomial-specific parameters
      list(tau_lower = binomial_tau_lower)
    }
  )
  stan_model <- get_stan_model(
    type = "glm",
    glm_family_string = cat_init$family,
    joint_tau = TRUE,
    glm_binomial_joint_theta = binomial_joint_theta,
    glm_binomial_joint_alpha = binomial_joint_alpha,
    glm_gaussian_known_variance = (!is.null(cat_init$custom_variance))
  )

  # Perform Bayesian Sampling via RStan
  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    ...
  )

  ## Extract the mean of coefs from the stan sampling model
  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]

  if (!is.null(coefs)) {
    names(coefs) <- c("(Intercept)", colnames(cat_init$adj_x))
  }

  if ((cat_init$family == "binomial") & (nrow(stan_sample_model) != 0)) {
    print_glm_bayes_joint_binomial_suggestion(
      alpha = alpha,
      stan_iter = iter,
      stan_sample_model = stan_sample_model,
      binomial_use_theta = binomial_use_theta,
      binomial_joint_alpha = binomial_joint_alpha,
      binomial_tau_lower = binomial_tau_lower,
    )
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm_bayes_joint",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau_alpha = ifelse(
      (binomial_joint_alpha & binomial_joint_theta),
      rstan::summary(stan_sample_model)$summary[
        grep("alpha", rownames(rstan::summary(stan_sample_model)$summary)),
        "mean"
      ],
      tau_alpha
    ),
    tau_gamma = tau_gamma,
    binomial_tau_lower = binomial_tau_lower,
    binomial_joint_theta = binomial_joint_theta,
    binomial_joint_alpha = binomial_joint_alpha,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    ...,
    ## Result
    model = stan_model,
    sample_model = stan_sample_model,
    coefficients = coefs,
    tau = if (binomial_joint_theta) {
      1 / rstan::summary(stan_sample_model)$summary[
        "theta", "mean"
      ]
    } else {
      rstan::summary(stan_sample_model)$summary[
        "tau", "mean"
      ]
    }
  )

  class(cat_model) <- c(cat_model$class, "cat_glm_bayes_joint", "bayes")

  return(cat_model)
}

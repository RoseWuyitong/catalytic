#' Bayesian GLM with Gibbs and HMC Sampling
#'
#' Fits a Bayesian generalized linear model (GLM) using a combination of Gibbs sampling
#' and Hamiltonian Monte Carlo (HMC) for sampling tau and coefs.
#'
#' @param formula The model formula specifying the response variable and predictors.
#' @param family GLM family to use. Defaults inferred from `cat_init`).
#' @param cat_init Initialized data object from `catalytic::cat_glm_initialization`.
#' @param iter Total number of Gibbs sampling iterations. Defaults to 2000.
#' @param gibbs_warmup Number of warmup iterations in Gibbs sampling. Defaults to 1000.
#' @param ... Additional arguments passed to `stats:glm`, for generating `coefs_0`
#' @param coefs_iter Number of iterations in the HMC sampling for coefs. Defaults to 5.
#' @param tau_0 Initial value for tau (downweighting factor for synthetic data). Defaults to one forth of the number of predictors.
#' @param coefs_0 Initial coefs for the GLM model. Defaults to simple glm model within function.
#' @param tau_alpha Shape parameter for the tau_gamma distribution used in tau sampling. Defaults to 2.
#' @param tau_gamma Rate parameter for the tau_gamma distribution used in tau sampling. Defaults to 1.
#' @param refresh Logical flag indicating whether to print progress in Gibbs sampling. Defaults to TRUE.
#'
#' @return An object of class `cat_glm_bayes_joint_gibbs` containing:
#'   \item{function_name}{Function name (`cat_glm_bayes_joint_gibbs`).}
#'   \item{formula}{Formula used for modeling.}
#'   \item{family}{GLM family used for modeling.}
#'   \item{cat_init}{Initialized data object from `cat_glm_initialization`.}
#'   \item{iter}{Total number of Gibbs sampling iterations.}
#'   \item{gibbs_warmup}{Number of warmup iterations for Gibbs sampling.}
#'   \item{coefs_iter}{Number of iterations for updating coefs in the HMC step.}
#'   \item{tau_0}{Initial value for tau used in modeling.}
#'   \item{coefs_0}{Initial coefs used in HMC (if provided).}
#'   \item{tau_alpha}{Shape parameter for the tau_gamma distribution in Gibbs sampling.}
#'   \item{tau_gamma}{Rate parameter for the tau_gamma distribution in Gibbs sampling.}
#'   \item{sys_time}{Timestamp of the sampling process completion.}
#'   \item{tau_model}{Fitted GLM model object using adjusted data and tau.}
#'   \item{inform_df}{Summary statistics of the Gibbs sampling output.}
#'   \item{gibbs_iteration_log}{Matrix containing the Gibbs sampling results for tau and coefs.}
#'   \item{tau}{Mean of sampled value of tau.}
#'   \item{coefs}{Mean of sampled coefs from the posterior distribution.}
#'
#' @export
#'
#' @examples
#' @references
#' Radford M. Neal. "MCMC using Hamiltonian dynamics." Handbook of Markov Chain Monte Carlo (2010).
#'
#' @details
#' This function fits a Bayesian GLM using Gibbs sampling to update tau and Hamiltonian
#' Monte Carlo (HMC) for updating coefs. It supports logistic regression models
#' (`binomial` family) and is designed to handle synthetic data in a Bayesian framework.
#'
#' Gibbs sampling iteratively updates tau and HMC samples coefs by simulating
#' Hamiltonian dynamics. The function provides flexibility in specifying initial values,
#' sampling iterations, and control over the sampling algorithms used.
#'
#' @note
#' - This function assumes `cat_init` is created using `catalytic::cat_glm_initialization`.
#' - Ensure `iter` is greater than `gibbs_warmup` to ensure proper sampling.
#'
#' @seealso \code{\link{HMC}}, \code{\link{cat_glm_initialization}}
cat_glm_bayes_joint_gibbs <- function(formula,
                                      cat_init,
                                      iter = 1000,
                                      warmup = 500,
                                      coefs_iter = 5,
                                      tau_0 = NULL,
                                      tau_alpha = 2,
                                      tau_gamma = 1,
                                      refresh = TRUE) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    gibbs_iter = iter,
    gibbs_warmup = warmup,
    coefs_iter = coefs_iter,
    tau_0 = tau_0,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  tau <- ifelse(is.null(tau_0), ncol(cat_init$adj_x) / 4, tau_0)

  kappa_ <- get_glm_log_density(
    family_string = cat_init$family,
    X = cat_init$adj_syn_x,
    Y = cat_init$syn_y,
    coefs = coef(do.call(
      stats::glm,
      list(
        formula = full_formula,
        family = cat_init$family,
        data = cat_init$syn_data
      )
    )),
    weights = 1 / cat_init$syn_size
  )

  # Fit GLM to obtain initial coefficients
  coefs <- coef(do.call(
    stats::glm,
    list(
      formula = full_formula,
      family = cat_init$family,
      data = cat_init$data,
      weights = c(
        rep(1, cat_init$obs_size),
        rep(tau / cat_init$syn_size, cat_init$syn_size)
      )
    )
  ))

  # Initialize a matrix to collect samples
  gibbs_iteration_log <- matrix(ncol = (1 + length(coefs)))
  colnames(gibbs_iteration_log) <- c("(Intercept)", colnames(cat_init$adj_x), "tau")

  if (refresh) {
    cat(
      "\nSAMPLING FROM GIBB'S NOW",
      "\nGibb's sampling:",
      "\nGibb's sampling: \n"
    )
  }

  # Define print boundaries for progress reporting
  print_boundary <- ifelse(
    iter >= 10,
    c(1, round(seq(0.1, 1, 0.1) * iter)),
    seq(1, iter)
  )

  # Initialize flags and timers for tracking progress
  sampling_start_flag <- FALSE
  warmup_start_time <- Sys.time()

  for (gibbs_i in 1:iter) {
    # Update tau by using the gamma distribution
    tau <- rgamma(
      n = 1,
      shape = ncol(cat_init$adj_x) + tau_alpha,
      scale = kappa_ + 1 / tau_gamma - get_glm_log_density(
        family = family,
        X = cat_init$adj_syn_x,
        Y = cat_init$syn_y,
        coefs = coefs,
        weights = 1 / cat_init$syn_size
      )
    )

    tau_model <- stats::glm(
      formula = full_formula,
      family = family,
      data = cat_init$data,
      weights = c(
        rep(1, cat_init$obs_size),
        rep(tau / cat_init$syn_size, cat_init$syn_size)
      )
    )

    # Compute diagonal approximation for HMC
    hmc_scale <- sqrt(get_glm_diag_approx_cov(cat_init$adj_x, tau_model))

    tau_weights <- c(
      rep(1, cat_init$obs_size),
      rep(tau / cat_init$syn_size, cat_init$syn_size)
    )

    # Update coefs by performing HMC sampling
    coefs <- get_hmc_mcmc_result(
      neg_log_den_func = function(coefs) {
        return(
          -get_glm_log_density(
            family = cat_init$family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefs,
            weights = tau_weights
          )
        )
      },
      neg_log_den_grad_func = function(coefs) {
        return(
          -get_glm_log_density_grad(
            family = cat_init$family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefs,
            weights = tau_weights
          )
        )
      },
      coefs_0 = coefs,
      iter = coefs_iter,
      hmc_scale = hmc_scale
    )

    gibbs_iteration_log[gibbs_i, ] <- c(coefs, tau)

    # Print progress updates
    if (((gibbs_i %in% print_boundary) | gibbs_i == warmup) & refresh) {
      print_str <- paste("Gibb's sampling:", gibbs_i, " / ", iter, "[", round(gibbs_i / iter * 100), "%]")

      if (gibbs_i <= warmup) {
        cat(print_str, " (Warmup)\n")
      } else {
        if (!sampling_start_flag) {
          sampling_start_time <- Sys.time()
          sampling_start_flag <- TRUE
        }
        cat(print_str, " (Sampling)\n")
      }
    }
  }

  sampling_end_time <- Sys.time()

  # Print elapsed times for warmup and sampling
  if (refresh) {
    cat(paste0(
      "Gibb's sampling: \n",
      "Gibb's sampling:",
      "  Elapsed Time: ",
      round(difftime(
        sampling_start_time,
        warmup_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Warm-up)\n",
      "Gibb's sampling:                ",
      round(difftime(
        sampling_end_time,
        sampling_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Sampling)\n",
      "Gibb's sampling:                ",
      round(difftime(
        sampling_end_time,
        warmup_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Total)\n",
      "Gibb's sampling: \n"
    ))

    # Finalize Setup and Output
    cat_model <- list(
      function_name = "cat_glm_bayes_joint_gibbs",
      ## Input/Processed parameters
      formula = formula,
      cat_init = cat_init,
      iter = iter,
      warmup = warmup,
      coefs_iter = coefs_iter,
      tau_0 = tau_0,
      tau_alpha = tau_alpha,
      tau_gamma = tau_gamma,
      refresh = refresh,
      sys_time = sampling_end_time,
      ## Result
      gibbs_iteration_log = gibbs_iteration_log,
      inform_df = tryCatch(
        {
          t(apply(
            gibbs_iteration_log,
            2,
            function(x) {
              c(
                mean = mean(x),
                se_mean = sd(x) / sqrt(length(x)),
                sd = sd(x),
                quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                n_eff = as.double(coda::effectiveSize(coda::as.mcmc(x)))
              )
            }
          ))
        },
        error = function(err) {
          t(sapply(
            gibbs_iteration_log,
            function(x) {
              c(
                mean = mean(x), se_mean = NA, sd = NA,
                quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), n_eff = 0
              )
            }
          ))
        }
      )
    )

    cat_model$tau <- cat_model$inform_df[(length(coefs) + 1), "mean"]
    cat_model$coefficients <- cat_model$inform_df[1:length(coefs), "mean"]

    class(cat_model) <- c(cat_model$class, "cat_glm_bayes_joint_gibbs", "bayes")
  }
  return(cat_model)
}

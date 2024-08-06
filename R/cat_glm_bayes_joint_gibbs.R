#' Bayesian GLM with Gibbs and HMC Sampling
#'
#' Fits a Bayesian generalized linear model (GLM) using a combination of Gibbs sampling
#' and Hamiltonian Monte Carlo (HMC) for sampling tau and coefficients.
#'
#' @param formula The model formula specifying the response variable and predictors.
#' @param family GLM family to use. Defaults inferred from `cat_init`).
#' @param cat_init Initialized data object from `catalytic::cat_glm_initialization`.
#' @param gibbs_iter Total number of Gibbs sampling iterations. Defaults to 2000.
#' @param gibbs_warmup Number of warmup iterations in Gibbs sampling. Defaults to 1000.
#' @param ... Additional arguments passed to `stats:glm`, for generating `init_coefficients`
#' @param coefficients_iter Number of iterations in the HMC sampling for coefficients. Defaults to 5.
#' @param init_tau Initial value for tau (downweighting factor for synthetic data). Defaults to one forth of the number of predictors.
#' @param init_coefficients Initial coefficients for the GLM model. Defaults to simple glm model within function.
#' @param alpha Shape parameter for the gamma distribution used in tau sampling. Defaults to 2.
#' @param gamma Rate parameter for the gamma distribution used in tau sampling. Defaults to 1.
#' @param stan_refresh Refresh rate for the Stan sampling. Defaults to 0.
#' @param gibbs_refresh Logical flag indicating whether to print progress in Gibbs sampling. Defaults to TRUE.
#'
#' @return An object of class `cat_glm_bayes_joint_gibbs` containing:
#'   \item{function_name}{Function name (`cat_glm_bayes_joint_gibbs`).}
#'   \item{formula}{Formula used for modeling.}
#'   \item{family}{GLM family used for modeling.}
#'   \item{cat_init}{Initialized data object from `cat_glm_initialization`.}
#'   \item{gibbs_iter}{Total number of Gibbs sampling iterations.}
#'   \item{gibbs_warmup}{Number of warmup iterations for Gibbs sampling.}
#'   \item{coefficients_iter}{Number of iterations for updating coefficients in the HMC step.}
#'   \item{init_tau}{Initial value for tau used in modeling.}
#'   \item{init_coefficients}{Initial coefficients used in HMC (if provided).}
#'   \item{alpha}{Shape parameter for the gamma distribution in Gibbs sampling.}
#'   \item{gamma}{Rate parameter for the gamma distribution in Gibbs sampling.}
#'   \item{sys_time}{Timestamp of the sampling process completion.}
#'   \item{coefs_tau_model}{Fitted GLM model object using adjusted data and tau.}
#'   \item{inform_df}{Summary statistics of the Gibbs sampling output.}
#'   \item{gibbs_collect}{Matrix containing the Gibbs sampling results for tau and coefficients.}
#'   \item{tau}{Mean of sampled value of tau.}
#'   \item{coefficients}{Mean of sampled coefficients from the posterior distribution.}
#'
#' @export
#'
#' @examples
#' @references
#' Radford M. Neal. "MCMC using Hamiltonian dynamics." Handbook of Markov Chain Monte Carlo (2010).
#'
#' @details
#' This function fits a Bayesian GLM using Gibbs sampling to update tau and Hamiltonian
#' Monte Carlo (HMC) for updating coefficients. It supports logistic regression models
#' (`binomial` family) and is designed to handle synthetic data in a Bayesian framework.
#'
#' Gibbs sampling iteratively updates tau and HMC samples coefficients by simulating
#' Hamiltonian dynamics. The function provides flexibility in specifying initial values,
#' sampling iterations, and control over the sampling algorithms used.
#'
#' @note
#' - This function assumes `cat_init` is created using `catalytic::cat_glm_initialization`.
#' - Ensure `gibbs_iter` is greater than `gibbs_warmup` to ensure proper sampling.
#'
#' @seealso \code{\link{HMC}}, \code{\link{cat_glm_initialization}}
cat_glm_bayes_joint_gibbs <- function(formula,
                                      family = NULL,
                                      cat_init,
                                      gibbs_iter = 1000,
                                      gibbs_warmup = 500,
                                      ...,
                                      coefficients_iter = 5,
                                      init_tau = NULL,
                                      init_coefficients = NULL,
                                      alpha = 2,
                                      gamma = 1,
                                      stan_refresh = 0,
                                      gibbs_refresh = TRUE) {
  if (!is.cat_glm_init(cat_init)) {
    stop(
      "Please enter `cat_init` that generated from `catalytic::cat_glm_initialization`."
    )
  }

  if (gibbs_iter <= gibbs_warmup) {
    stop(
      "Please note that gibbs_iter should be larger than gibbs_warmup."
    )
  }

  if (coefficients_iter > 10) {
    response <- menu(c("YES", "NO"), title = paste0(
      "Please note that large stan iteration times can",
      "slow down gibb's sampling process.",
      "\nDo you want to continue?\n"
    ))

    if (response == 2) {
      stop("Process terminated by user.", call. = FALSE)
    }
  }


  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  f_lhs <- get_formula_lhs(formula)
  if (f_lhs != "NULL" & f_lhs != colnames(cat_init$y)) stop("Please enter the correct y name for formula.")

  cat_init <- get_adjusted_data(cat_init = cat_init, formula = formula)

  if (is.null(family)) family <- cat_init$family
  tau <- if (is.null(init_tau)) ncol(cat_init$adj_x) / 4 else init_tau
  if (tau < 0) stop("Please enter non-negative `tau` value.")

  full_formula <- as.formula(paste0(colnames(cat_init$y), f_rhs))

  if (get_glm_family(family) != "binomial") {
    stop(
      "Please note that only logistic regression can use `cat_glm_bayes_joint_gibbs`."
    )
  }

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

  tau <- if (is.null(init_tau)) ncol(cat_init$adj_x) / 4 else init_tau

  suppressWarnings(coefs_tau_model <- do.call(stats::glm, list(
    formula = full_formula,
    family = family,
    data = cat_init$data,
    weights = c(
      rep(1, cat_init$obs_size),
      rep(tau / cat_init$syn_size, cat_init$syn_size)
    )
  )))

  hmc_scale <- sqrt(get_glm_diag_approx_cov(cat_init$adj_x, coefs_tau_model))

  coefficients <- if (is.null(init_coefficients)) coefs_tau_model$coefficients else init_coefficients

  gibbs_collect <- matrix(ncol = (1 + length(coefficients)))
  colnames(gibbs_collect) <- c("(Intercept)", colnames(cat_init$x), "tau")

  if (gibbs_refresh) {
    cat(
      "\nSAMPLING FROM GIBB'S NOW",
      "\nGibb's sampling:",
      "\nGibb's sampling: \n"
    )
  }


  print_boundary <- if (gibbs_iter >= 10) {
    c(
      1,
      round(seq(0.1, 1, 0.1) *
        gibbs_iter)
    )
  } else {
    seq(1, gibbs_iter)
  }


  sampling_start_flag <- TRUE
  warmup_start_time <- Sys.time()

  for (gibbs in 1:gibbs_iter) {
    # Update tau

    tau <- rgamma(
      n = 1,
      shape = ncol(cat_init$adj_x) + alpha,
      scale = kappa_ + 1 / gamma - get_glm_log_density(
        family = family,
        X = cat_init$adj_syn_x,
        Y = cat_init$syn_y,
        coefs = coefficients,
        weights = 1 / cat_init$syn_size
      )
    )


    # Update coefficients
    suppressWarnings(coefs_tau_model <- do.call(stats::glm, list(
      formula = full_formula,
      family = family,
      data = cat_init$data,
      weights = c(
        rep(1, cat_init$obs_size),
        rep(tau / cat_init$syn_size, cat_init$syn_size)
      ),
      ...
    )))

    hmc_scale <- sqrt(get_glm_diag_approx_cov(cat_init$adj_x, coefs_tau_model))

    weights <- c(
      rep(1, cat_init$obs_size),
      rep(tau / cat_init$syn_size, cat_init$syn_size)
    )

    coefficients <- hmc_mcmc(
      neg_log_den_func = function(coefficients) {
        return(
          -get_glm_log_density(
            family = family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefficients,
            weights = weights
          )
        )
      },
      neg_log_den_grad_func = function(coefficients) {
        return(
          -get_glm_log_density_grad(
            family = family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefficients,
            weights = weights
          )
        )
      },
      init_coefficients = coefficients,
      iter = coefficients_iter,
      hmc_scale = hmc_scale
    )



    if (gibbs > gibbs_warmup) {
      gibbs_collect <- rbind(gibbs_collect, c(coefficients, tau))

      if (sampling_start_flag) {
        sampling_start_time <- Sys.time()
        sampling_start_flag <- FALSE
      }

      if (gibbs %in% print_boundary & gibbs_refresh) {
        cat(
          "Gibb's sampling: ",
          gibbs,
          "/",
          gibbs_iter,
          " [",
          round(gibbs / gibbs_iter * 100),
          "%] ",
          " (Sampling)\n"
        )
      }
    } else if (((gibbs < gibbs_warmup & (gibbs %in% print_boundary)) | gibbs == gibbs_warmup) & gibbs_refresh) {
      cat(
        "Gibb's sampling: ",
        gibbs,
        "/",
        gibbs_iter,
        " [",
        round(gibbs / gibbs_iter * 100),
        "%] ",
        " (Warmup)\n"
      )
    }
  }

  sampling_end_time <- Sys.time()
  gibbs_collect <- gibbs_collect[-1, ]

  cat_glm_bayes_joint_gibbs_ <- list(
    function_name = "cat_glm_bayes_joint_gibbs",
    formula = formula,
    family = family,
    cat_init = cat_init,
    gibbs_iter = gibbs_iter,
    gibbs_warmup = gibbs_warmup,
    coefficients_iter = coefficients_iter,
    init_tau = init_tau,
    init_coefficients = init_coefficients,
    alpha = alpha,
    gamma = gamma,
    sys_time = sampling_end_time,
    coefs_tau_model = coefs_tau_model,
    gibbs_collect = gibbs_collect,
    inform_df = tryCatch(
      {
        t(apply(
          gibbs_collect,
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
          gibbs_collect,
          function(x) {
            c(
              mean = mean(x), se_mean = NA, sd = NA,
              quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), n_eff = 0
            )
          }
        ))
      }
    ),
    ...
  )

  cat_glm_bayes_joint_gibbs_$tau <- cat_glm_bayes_joint_gibbs_$inform_df[(length(coefficients) + 1), "mean"]
  cat_glm_bayes_joint_gibbs_$coefficients <- cat_glm_bayes_joint_gibbs_$inform_df[1:length(coefficients), "mean"]

  class(cat_glm_bayes_joint_gibbs_) <- c(cat_glm_bayes_joint_gibbs_$class, "cat_glm_bayes_joint_gibbs", "bayes")

  cat(paste0(
    "Gibb's sampling: \n",
    "Gibb's sampling:",
    "  Elapsed Time: ",
    round(difftime(sampling_start_time,
      warmup_start_time,
      units = "secs"
    ), 3),
    " seconds ",
    "(Warm-up)\n",
    "Gibb's sampling:                ",
    round(difftime(sampling_end_time,
      sampling_start_time,
      units = "secs"
    ), 3),
    " seconds ",
    "(Sampling)\n",
    "Gibb's sampling:                ",
    round(difftime(sampling_end_time,
      warmup_start_time,
      units = "secs"
    ), 3),
    " seconds ",
    "(Total)\n",
    "Gibb's sampling: \n"
  ))

  return(cat_glm_bayes_joint_gibbs_)
}


hmc_mcmc <- function(neg_log_den_func,
                     neg_log_den_grad_func,
                     init_coefficients,
                     iter = 5,
                     hmc_scale = 0.01,
                     hmc_steps = 5,
                     ...,
                     return_only_last = TRUE) {
  mc_collect <- matrix(0, iter + 1, length(init_coefficients))
  mc_collect[1, ] <- init_coefficients

  for (i in 1:iter) {
    mc_collect[i + 1, ] <- HMC(
      U = neg_log_den_func,
      grad_U = neg_log_den_grad_func,
      # add random value to epsilon to overcome poor performance
      epsilon = rbeta(1, 0.5, 0.5) * hmc_scale / hmc_steps,
      L = hmc_steps,
      current_q = mc_collect[i, ]
    )[[1]]
  }

  if (return_only_last) {
    return(mc_collect[iter + 1, ])
  }
  return(mc_collect[-1, ])
}


# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
#
# Radford M. Neal, 2010.
#
# This program appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# to appear in the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   U          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   grad_U     A function to evaluate the gradient of U.
#
#   epsilon    The stepsize to use for the leapfrog steps.
#
#   L          The number of leapfrog steps to do to propose a new state.
#
#   current_q  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_q if the endpoint of the
# trajectory was rejected).
#
# This function was written for illustrative purposes.  More elaborate
# implementations of this basic HMC method and various variants of HMC
# are available from my web page, http://www.cs.utoronto.ca/~radford/
HMC <- function(U, grad_U, epsilon, L, current_q) {
  q <- current_q
  p <- rnorm(length(q), 0, 1) # independent standard normal variates
  current_p <- p

  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2

  # Alternate full steps for position and momentum

  for (i in 1:L)
  {
    # Make a full step for the position

    q <- q + epsilon * p

    # Make a full step for the momentum, except at end of trajectory

    if (i != L) p <- p - epsilon * grad_U(q)
  }

  # Make a half step for momentum at the end.

  p <- p - epsilon * grad_U(q) / 2

  # Negate momentum at end of trajectory to make the proposal symmetric

  p <- -p

  # Evaluate potential and kinetic energies at start and end of trajectory

  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
    return(list(q, proposed_U, accept = 1)) # accept
  } else {
    return(list(current_q, current_U, accept = 0)) # reject
  }
}

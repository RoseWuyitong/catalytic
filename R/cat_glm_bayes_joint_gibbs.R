#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param ...
#' @param init_tau
#' @param init_beta
#' @param gibbs_iter
#' @param algorithm
#' @param alpha
#' @param gamma
#' @param gibbs_warmup
#' @param beta_iter
#' @param gibbs_print
#' @param stan_refresh
#'
#' @return
#' @export
#'
#' @examples
cat_glm_bayes_joint_gibbs <- function(formula,
                                      family = binomial(),
                                      data_,
                                      init_tau = NULL,
                                      init_coefficients = NULL,
                                      gibbs_iter = 1000,
                                      gibbs_warmup = 500,
                                      ...,
                                      coefficients_iter = 5,
                                      algorithm = "HMC",
                                      alpha = 2,
                                      gamma = 1,
                                      stan_refresh = 0,
                                      gibbs_print = TRUE) {
  .not_true_stop(
    is.data_(data_),
    paste(
      "Please enter `data_` that generated from",
      "catalytic::data_initialization."
    )
  )

  .small_number_stop(
    gibbs_iter,
    gibbs_warmup + 1,
    "Please note that gibbs_iter should be larger than gibbs_warmup."
  )

  .large_number_ask(
    coefficients_iter,
    10,
    paste(
      "Please note that large stan iteration times can",
      "slow down gibb's sampling process, are you sure?\ny/n"
    ),
    "n"
  )

  data_ <- get_adj_data(data_ = data_, formula = formula)

  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- get_family(family)

  .not_same_stop(
    fam,
    "binomial",
    paste(
      "Please note that only logistic regression can use,",
      "`cat_glm_bayes_joint_gibbs`."
    )
  )

  print_boundary <- if (gibbs_iter >= 10) {
    c(
      1,
      round(seq(0.1, 1, 0.1) *
        gibbs_iter)
    )
  } else {
    seq(1, gibbs_iter)
  }

  kappa_ <- get_log_density(fam)(
    X = data_$adj_syn_x,
    Y = data_$syn_y,
    model = do.call(stats::glm, list(formula = formula, data = data_$syn_data, family = family)),
    weights = 1 / data_$syn_size)

  tau <- if (is.null(init_tau)) ncol(data_$adj_x) / 4 else init_tau

  .small_number_stop(tau, 0, "Please enter non-negative `init_tau` value.")

  gamma_function <- (ncol(data_$adj_x) + alpha - 1) * log(tau) -
    (tau) * (kappa_ + gamma^-1)

  sampling_start_time_flag <- TRUE
  warmup_start_time <- Sys.time()
  tau <- tau
  sample_model <- do.call(stats::glm, list(
    formula = formula,
    family = family,
    data = data_$data,
    weights = c(
      rep(1, data_$obs_size),
      rep(tau / data_$syn_size, data_$syn_size)
    )
  ))

  hmc_scale <- sqrt(get_diag_approx_cov(data_$adj_x, sample_model))

  coefficients <- if (is.null(init_coefficients)) sample_model$coefficients else init_coefficients

  gibbs_collect <- matrix(ncol = (1 + length(coefficients)))
  colnames(gibbs_collect) <- c("tau", paste0("beta[", 1:length(coefficients), "]"))

  if (gibbs_print) {
    cat(
      "\nSAMPLING FROM GIBB'S NOW",
      "\nGibb's sampling:",
      "\nGibb's sampling: \n"
    )
  }

  for (gibbs in 1:gibbs_iter) {
    # Update tau
    tau <- rgamma(
      n = 1,
      shape = ncol(data_$adj_x) + alpha,
      scale = kappa_ +
        1 / gamma -
        get_log_density(fam)(
          data_$adj_syn_x,
          data_$syn_y,
          coefficients = coefficients,
          weights = 1 / data_$syn_size)
    )

    # Update coefficients
    sample_model <- do.call(stats::glm, list(
      formula = formula,
      family = family,
      data = data_$data,
      weights = c(
        rep(1, data_$obs_size),
        rep(tau / data_$syn_size, data_$syn_size)
      )
    ))

    hmc_scale <- sqrt(get_diag_approx_cov(data_$adj_x, sample_model))

    weights <- c(
      rep(1, data_$obs_size),
      rep(tau / data_$syn_size, data_$syn_size)
    )

    coefficients <- hmc_mcmc(
      neg_log_den_func = function(coefficients) {
        return(-get_log_density(fam)(
          X = data_$adj_x,
          Y = data_$y,
          coefficients = coefficients,
          weights = weights
        ))
      },
      neg_log_den_grad_func = function(coefficients) {
        return(-get_log_density_grad(fam)(
          X = data_$adj_x,
          Y = data_$y,
          coefficients = coefficients,
          weights = weights
        ))
      },
      init_coefficients = coefficients,
      iter = coefficients_iter,
      hmc_scale = hmc_scale
    )



    if (gibbs > gibbs_warmup) {
      gibbs_collect <- rbind(gibbs_collect, c(tau, coefficients))

      if (sampling_start_time_flag) {
        sampling_start_time <- Sys.time()
        sampling_start_time_flag <- FALSE
      }

      if (gibbs %in% print_boundary & gibbs_print) {
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
    } else if (((gibbs < gibbs_warmup & (gibbs %in% print_boundary)) | gibbs == gibbs_warmup) & gibbs_print) {
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

  cat_gibbs <- list(
    cat_function = "cat_glm_bayes_joint_gibbs",
    formula = formula,
    data_ = data_,
    family = family,
    algorithm = algorithm,
    sample_model = sample_model,
    gibbs_collect = gibbs_collect,
    gibbs_iter = gibbs_iter,
    gibbs_warmup = gibbs_warmup,
    alpha = alpha,
    gamma = gamma,
    coefficients_iter = coefficients_iter,
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
    )
  )

  cat_gibbs$tau <- cat_gibbs$inform_df[1, "mean"]
  cat_gibbs$coefficients <- cat_gibbs$inform_df[-1, "mean"]

  class(cat_gibbs) <- c(cat_gibbs$class, "bayes_gibbs", "cat_glm_bayes_joint_gibbs")

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
add_gibbs_iter.cat_gibbs_glm <- function(cat_gibbs, gibbs_iter = 100) {
  updated_cat_gibbs <- cat_gibbs_glm(
    formula = cat_gibbs$formula,
    family = cat_gibbs$family,
    data_ = cat_gibbs$data_,
    init_tau = cat_gibbs$gibbs_collect[nrow(cat_gibbs$gibbs_collect), "tau"],
    init_coefficients = cat_gibbs$gibbs_collect[
      nrow(cat_gibbs$gibbs_collect),
      grep("beta", colnames(cat_gibbs$gibbs_collect))
    ],
    gibbs_iter = gibbs_iter,
    gibbs_warmup = 0,
    algorithm = cat_gibbs$algorithm,
    alpha = cat_gibbs$alpha,
    gamma = cat_gibbs$gamma,
    coefficients_iter = cat_gibbs$coefficients_iter
  )

  updated_cat_gibbs$gibbs_iter <- cat_gibbs$gibbs_iter + gibbs_iter
  updated_cat_gibbs$gibbs_collect <- rbind(cat_gibbs$gibbs_collect, updated_cat_gibbs$gibbs_collect)

  updated_cat_gibbs$inform_df <- t(apply(
    updated_cat_gibbs$gibbs_collect,
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
  updated_cat_gibbs$tau <- updated_cat_gibbs$inform_df[1, "mean"]
  updated_cat_gibbs$coefficients <- updated_cat_gibbs$inform_df[-1, "mean"]
  return(updated_cat_gibbs)
}

add_gibbs_iter <- function(cat_gibbs, gibbs_iter) {
  UseMethod("add_gibbs_iter")
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

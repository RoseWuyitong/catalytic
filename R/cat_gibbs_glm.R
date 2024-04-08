#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param ...
#' @param iter
#' @param init_tau
#' @param init_beta
#' @param gibbs_iter
#' @param algorithm
#' @param alpha
#' @param gamma
#' @param stan_refresh
#' @param warmup
#'
#' @return
#' @export
#'
#' @examples
cat_gibbs_glm <- function(formula,
                          family = gaussian(),
                          data_,
                          init_tau = NULL,
                          init_beta = NULL,
                          gibbs_iter = 2000,
                          gibbs_warmup = 1000,
                          ...,
                          stan_iter = 3,
                          stan_warmup = 2,
                          algorithm = "HMC",
                          alpha = 2,
                          gamma = 1,
                          stan_refresh = 0,
                          gibbs_print = TRUE) {
  .not_in_stop(
    c(
      "obs_size", "obs_data", "obs_x", "obs_y", "syn_size", "syn_data", "syn_x",
      "syn_y", "syn_model", "data", "size", "x", "y"
    ),
    names(data_),
    "Not all necessary variables are present in data_."
  )

  .small_number_stop(
    gibbs_iter,
    gibbs_warmup + 1,
    "gibbs_iter should be larger than gibbs_warmup."
  )

  .large_number_ask(
    stan_iter,
    10,
    "large stan iteration times can slow down gibb's sampling process,
  are you sure? y/n",
    "n"
  )

  data_$adj_x <- model.frame(rhs(formula, TRUE), data_$x)
  data_$adj_obs_x <- data_$adj_x[1:data_$obs_size, ]
  data_$adj_syn_x <- data_$adj_x[(data_$obs_size + 1):data_$size, ]

  validate_glm_formula(formula)
  family <- validate_family(family)
  fam <- get_family(family)

  print_boundary <- if (gibbs_iter >= 10) c(1, round(seq(0.1, 1, 0.1) * gibbs_iter)) else seq(1, gibbs_iter)

  kappa_ <- log_density(fam)(
    data_$adj_syn_x,
    data_$syn_y,
    stats::glm(formula = formula, data = data_$syn_data, family = family),
    weights = 1 / data_$syn_size)

  data <- list(
    obs_n = data_$obs_size,
    syn_n = data_$syn_size,
    x_col_n = ncol(data_$adj_x),
    tau = if (is.null(init_tau)) ncol(data_$adj_x) / 4 else init_tau,
    obs_x = data_$adj_obs_x,
    obs_y = as.numeric(data_$obs_y[[1]]),
    syn_x = data_$adj_syn_x,
    syn_y = as.numeric(data_$syn_y[[1]]),
    kappa = kappa_
  )

  stan_model <- stan_model_selection(
    fam = fam,
    joint = TRUE,
    gibbs = TRUE,
    use_theta = FALSE
  )


  sampling_start_time_flag <- TRUE
  warmup_start_time <- Sys.time()
  if (is.null(init_beta)) {
    suppressWarnings(sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      chains = 1,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm,
      refresh = stan_refresh
    ))

    init_beta <- rstan::extract(sample_model)$beta[
      nrow(rstan::extract(sample_model)$beta),
    ]
  }

  gibbs_collect <- matrix(ncol = (1 + length(init_beta)))
  colnames(gibbs_collect) <- c("tau", paste0("beta[", 1:length(init_beta), "]"))

  if (gibbs_print) {
    cat(
      "\nSAMPLING FROM GIBB'S NOW",
      "\nGibb's sampling:",
      "\nGibb's sampling: \n") }

  for (gibbs in 1:gibbs_iter) {
    suppressWarnings(sample_model <- rstan::sampling(
      object = stan_model,
      data = data,
      init = function(...) list(beta = init_beta),
      chains = 1,
      iter = iter,
      warmup = warmup,
      algorithm = algorithm,
      refresh = stan_refresh
    ))

    init_beta <- rstan::extract(sample_model)$beta[
      nrow(rstan::extract(sample_model)$beta),
    ]

    data$tau <- rgamma(
      n = 1,
      shape = ncol(data_$adj_x) + alpha,
      scale = kappa_ +
        1 / gamma -
        log_density(fam)(
          data_$adj_syn_x,
          data_$syn_y,
          beta = init_beta,
          weights = 1 / data_$syn_size)
    )


    if (gibbs > gibbs_warmup) {
      gibbs_collect <- rbind(gibbs_collect, c(data$tau, init_beta))

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
    cat_function = "cat_gibbs_glm",
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
    stan_iter = stan_iter,
    stan_warmup = stan_warmup,
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
  cat_gibbs$beta <- cat_gibbs$inform_df[-1, "mean"]

  class(cat_gibbs) <- c(cat_gibbs$class, "cat_gibbs_glm", "cat_glm")

  cat(paste0(
    "Gibb's sampling: \n",
    "Gibb's sampling:", "  Elapsed Time: ", round(sampling_start_time - warmup_start_time, 3), " seconds ", "(Warm-up)\n",
    "Gibb's sampling:                ", round(sampling_end_time - sampling_start_time, 3), " seconds ", "(Sampling)\n",
    "Gibb's sampling:                ", round(sampling_end_time - warmup_start_time, 3), " seconds ", "(Total)\n",
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
    init_beta = cat_gibbs$gibbs_collect[
      nrow(cat_gibbs$gibbs_collect),
      grep("beta", colnames(cat_gibbs$gibbs_collect))
    ],
    gibbs_iter = gibbs_iter,
    gibbs_warmup = 0,
    algorithm = cat_gibbs$algorithm,
    alpha = cat_gibbs$alpha,
    gamma = cat_gibbs$gamma,
    stan_iter = cat_gibbs$stan_iter,
    stan_warmup = cat_gibbs$stan_warmup
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
          )}))
  updated_cat_gibbs$tau <- updated_cat_gibbs$inform_df[1, "mean"]
  updated_cat_gibbs$beta <- updated_cat_gibbs$inform_df[-1, "mean"]
  return(updated_cat_gibbs)
}

add_gibbs_iter <- function(cat_gibbs, gibbs_iter) {
  UseMethod("add_gibbs_iter")
}

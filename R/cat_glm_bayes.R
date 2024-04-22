#' Title
#'
#' @param formula
#' @param family
#' @param data_
#' @param tau
#' @param ...
#' @param chains
#' @param iter
#' @param algorithm
#' @param warmup
#'
#' @return
#' @export
#'
#' @examples
cat_glm_bayes <- function(formula,
                          family = gaussian(),
                          data_,
                          tau = NULL,
                          ...,
                          chains = 4,
                          iter = 2000,
                          warmup = 1000,
                          algorithm = "HMC",
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

  kappa_ <- log_density(fam)(
    X = data_$adj_syn_x,
    Y = data_$syn_y,
    model = stats::glm(formula = formula, data = data_$syn_data, family = family),
    weights = 1 / data_$syn_size)

  data <- list(
    obs_n = data_$obs_size,
    syn_n = data_$syn_size,
    x_col_n = ncol(data_$adj_x),
    tau = if (is.null(tau)) ncol(data_$adj_x) / 4 else tau,
    obs_x = data_$adj_obs_x,
    obs_y = as.numeric(data_$obs_y[[1]]),
    syn_x = data_$adj_syn_x,
    syn_y = as.numeric(data_$syn_y[[1]])
  )

  .small_number_stop(data$tau, 0, "Please enter non-negative `tau` value.")

  stan_model <- stan_model_selection(
    fam = fam,
    joint = FALSE,
    gibbs = FALSE,
    use_theta = FALSE
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

  cat_glm_obj <- list(
    function_name = "cat_glm_bayes",
    formula = formula,
    data_ = data_,
    tau = data$tau,
    family = family,
    algorithm = algorithm,
    sample_model = sample_model,
    coefficients = rstan::summary(sample_model)$summary[
      grep("beta", rownames(rstan::summary(sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_glm_obj) <- c(cat_glm_obj$class, "bayes", "cat_glm_bayes")

  return(cat_glm_obj)
}

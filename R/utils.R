log_density <- function(family) {
  switch(family,
    "binomial" = function(X, Y, model, weights = FALSE, ...) {
      theta <- matrix_multiply(X, model$coefficients)
      if (isFALSE(weights)) {
        return(sum(theta * Y - log(1 + exp(theta))))
      } else {
        return(weights * sum(theta * Y - log(1 + exp(theta))))
      }
    },
    "gaussian" = function(X, Y, model, weights = FALSE, ...) {
      if (isFALSE(weights)) {
        return(0.5 * (1 + log(2 * pi * var(model$residuals))))
      } else {
        return(0.5 * weights * (1 + log(2 * pi * var(model$residuals))))
      }
    }
  )
}

fam_mean <- function(family) {
  switch(family,
    "binomial" = function(X, beta, ...) 1 / (1 + exp(-matrix_multiply(X, beta))),
    "gaussian" = function(X, beta, ...) matrix_multiply(X, beta),
    stop("please enter right regression type.")
  )
}

y_value <- function(family) {
  switch(family,
    "binomial" = function(X, model) {
      generate_data_from_dist(
        "binomial",
        nrow(X),
        fam_mean("binomial")(X, model$coefficients)
      )
    },
    "gaussian" = function(X, model) {
      generate_data_from_dist(
        "gaussian",
        nrow(X),
        fam_mean("gaussian")(X, model$coefficients),
        sqrt(sum(model$residuals^2) / model$df.residual)
      )
    }
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, family, ...) 2 * fam_mean(family)(X, beta),
    "classification_error" = function(X, beta, family, ...) 2 * fam_mean(family)(X, beta),
    "logistic_deviance" = function(X, beta, family = "binomial", ...) {
      if (!is.binomial(family)) {
        stop("only binomial family can use logistic_deviance", call. = FALSE)
      } else {
        matrix_multiply(X, beta)
      }
    }
  )
}

correlation_matrix <- function(corr, dim) matrix(corr, nrow = dim, ncol = dim) + diag(dim) * (1 - corr)

standard_error <- function(X) sd(X) / sqrt(length(X))

matrix_multiply <- function(X, Y, intl = TRUE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (intl) {
    if (ncol(X) + 1 == nrow(Y)) {
      return(Y[1] + X %*% Y[-1])
    } else if (
      ncol(Y) + 1 == nrow(X)) {
      return(X[1] + Y %*% X[-1])
    } else if (
      nrow(X) + 1 == nrow(Y)) {
      return(Y[1] + Y[-1] %*% X)
    } else if (
      nrow(Y) + 1 == nrow(X)) {
      return(X[1] + X[-1] %*% Y)
    }
  }

  if (ncol(X) == nrow(Y)) {
    return(X %*% Y)
  } else if (
    nrow(X) == ncol(Y)) {
    return(Y %*% X)
  }

  stop(paste("Cannot multiply these two matrix", str(X), str(Y)))
}

adjust_data_type <- function(X) {
  X <- as.data.frame(X)
  for (i in 1:ncol(X)) {
    if (!(NA %in% suppressWarnings(as.numeric(X[, i])))) {
      X[, i] <- as.numeric(X[, i])
    } else if (!(NA %in% suppressWarnings(as.logical(X[, i])))) {
      X[, i] <- as.logical(X[, i])
    } else {
      X[, i] <- as.factor(X[, i])
    }
  }
  return(X)
}

get_iqr <- function(x) {
  right_quantile <- 0.75
  iqr_value <- 0
  while (iqr_value == 0 & right_quantile < 1) {
    iqr_value <- diff(quantile(as.numeric(x),
      c(1 - right_quantile, right_quantile),
      na.rm = FALSE,
      names = FALSE,
      type = 7
    ))
    right_quantile <- right_quantile + 0.05
  }

  return(iqr_value)
}

generate_data_from_dist <- function(family = "gaussian",
                                    sample_size = 10,
                                    mean_ = 0,
                                    sd = 1) {
  switch(family,
    "binomial" = rbinom(
      sample_size,
      1,
      mean_
    ),
    "gaussian" = rnorm(
      sample_size,
      mean_,
      sd
    ),
    "uniform" = rep(mean_, sample_size)
  )
}

discrepancy_functions <- function(discrepancy_function) {
  switch(discrepancy_function,
    "logarithmic" = function(x, y, model, family) {
      # cap the value between 0.0001 and 0.9999 for better calculation process
      x_cap <- pmax(0.0001, pmin(0.9999, fam_mean(family)(x, model$coefficients)))
      y_cap <- pmax(0.0001, pmin(0.9999, y))
      x * log(x_cap / y_cap) + (1 - x) * log((1 - x_cap) / (1 - y_cap))
    },
    "square_error" = function(x, y, model, family) (y - fam_mean(family)(x, model$coefficients))^2,
    "classification_error" = function(x, y, model, family) as.numeric(y == fam_mean(family)(x, model$coefficients)),
    "logistic_deviance" = function(x, y, model, family) -log_density(family)(x, y, model)
  )
}


lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, family) 2 * fam_mean(family)(X, beta),
    "classification_error" = function(X, beta, family) 2 * fam_mean(family)(X, beta),
    "logistic_deviance" = function(X, beta, family) matrix_multiply(X, beta)
  )
}

gelman_rubin_statistic <- function(rstan_model, rhat_threhold = 0.01) {
  if (all(
    abs(
      head(
        rstan::summary(rstan_model)$summary[, "Rhat"], -1
      ) - 1
    ) < rhat_threhold
  )) {
    return("PASS")
  }

  warning("Burn-in//warm up times are not large enought.")
  return("FAIL")
}


effective_sample_size <- function(rstan_model, ess_threhold = 100) {
  if (all(
    head(
      rstan::summary(rstan_model)$summary[, "n_eff"], -1
    ) > ess_threhold
  )) {
    return("PASS")
  }

  warning("Iteration times are not large enought.")
  return("FAIL")
}

standardize_data <- function(data, if_missing = NULL) {
  if (missing(data) || is.null(data)) {
    warn_data_arg_missing()
    return(if_missing)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  drop_redundant_dims(data)

  for (i in 1:ncol(data)) {
    if (!(NA %in% suppressWarnings(as.numeric(data[, i])))) {
      data[, i] <- as.numeric(data[, i])
    } else if (!(NA %in% suppressWarnings(as.logical(data[, i])))) {
      data[, i] <- as.logical(data[, i])
    } else {
      data[, i] <- as.factor(data[, i])
    }
  }

  return(data)
}


drop_redundant_dims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop = FALSE], drop)
  return(data)
}

validate_glm_formula <- function(f) {
  if (!any(grepl("~", f))) {
    stop("please enter ~ for formula",
      call. = FALSE
    )
  }
}

validate_tau <- function(t) {
  if (!all(t >= 0)) {
    stop("please enter positive tau value",
      call. = FALSE
    )
  }
}



is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}
is.cat_obj <- function(x) inherits(x, "cat_obj")
is.initial_obj <- function(x) inherits(x, "initial_obj")

validate_family <- function(f) {
  if (is.character(f)) {
    f <- get(f, mode = "function", envir = parent.frame(2))
  }
  if (is.function(f)) {
    f <- f()
  }
  if (!is(f, "family")) {
    stop("'family' must be a family.", call. = FALSE)
  }

  return(f)
}

rhs <- function(formula, with_ = FALSE) {
  # Split the formula string by "~" and extract the first element
  if (with_) {
    return(as.formula(paste("~", deparse(rlang::f_rhs(formula)), collapse = " ")))
  }
  return(rlang::f_rhs(formula))
}

lhs <- function(formula) {
  return(rlang::f_lhs(formula))
}

.large_number_warning <- function(num, threshold, suggestion) {
  if (num > threshold) warning(suggestion, call. = FALSE)
}

.small_number_warning <- function(num, threshold, suggestion) {
  if (num < threshold) warning(suggestion, call. = FALSE)
}

.not_same_stop <- function(num1, num2, suggestion) {
  if (num1 != num2) stop(suggestion, call. = FALSE)
}

.not_in_stop <- function(obj, expected_objs, suggestion) {
  if (!all(obj %in% expected_objs)) stop(suggestion, call. = FALSE)
}

.not_single_str_stop <- function(obj, suggestion) {
  if (!is.character(obj) | length(obj) != 1) stop(suggestion, call. = FALSE)
}

.large_number_ask <- function(num1, threshold, question, stop_ans) {
  if (num1 > threshold) {
    user_input <- readline(prompt = question)
    if (identical(user_input, stop_ans)) stop("program stop", call. = FALSE)
  }
}


log_density <- function(family) {
  switch(family,
    "binomial" = function(X, Y, model, weights = FALSE, ...) {
      theta <- matrix_multiply(X, model$coefficients)
      if (isFALSE(weights)) {
        return(-sum(theta * Y - log(1 + exp(theta))))
      } else {
        return(-weights * sum(theta * Y - log(1 + exp(theta))))
      }
    },
    "gaussian" = function(X, Y, model, weights = FALSE, ...) {
      theta <- matrix_multiply(X, model$coefficients)
      if (isFALSE(weights)) {
        return(- 0.5 * (1+log(2*pi*var(model$residuals))))
      } else {
        return(- 0.5 * weights * (1+log(2*pi*var(model$residuals)))) # TODO: check
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
        sqrt(sum(model$residuals^2) / model$df.residual) # TODO: check
      )
    }
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, family, ...) 2 * fam_mean(family)(X, beta),
    "classification_error" = function(X, beta, family, ...) 2 * fam_mean(family)(X, beta),
    "logistic_deviance" = function(X, beta, family = "binomial", ...) matrix_multiply(X, beta)
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


fit <- function(x,
                y,
                weights = NULL,
                method = "mle",
                family = "gaussian") {
  switch(method,
    "mle" = if (identical(x, 1)) {
      stats::glm(y ~ 1, family = family, weights = weights)
    } else {
      stats::glm(y ~ ., data = as.data.frame(x), family = family, weights = weights)
    },
    "pseudo_mle" = if (identical(x, 1)) {
      stats::glm(y ~ 1, family = family, weights = weights, maxit = 25)
    } else {
      stats::glm(y ~ ., data = as.data.frame(x), family = family, weights = weights, maxit = 25)
    },
    "cauchy" = if (identical(x, 1)) {
      arm::bayesglm(y ~ 1, family = family, weights = weights)
    } else {
      arm::bayesglm(y ~ ., data = as.data.frame(x), family = family, weights = weights)
    }
  )
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

generate_data_from_dist <- function(generate_dist = "gaussian",
                                    sample_size = 10,
                                    mean_ = 0,
                                    sd = 1) {
  switch(generate_dist,
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
    "logarithmic" = function(x, y, beta, family) {
      # cap the value between 0.0001 and 0.9999 for better calculation process
      x_cap <- pmax(0.0001, pmin(0.9999, fam_mean(family)(x, beta)))
      y_cap <- pmax(0.0001, pmin(0.9999, y))
      x * log(x_cap / y_cap) + (1 - x) * log((1 - x_cap) / (1 - y_cap))
    },
    "square_error" = function(x, y, beta, family) (y - fam_mean(family)(x, beta))^2,
    "classification_error" = function(x, y, beta, family) as.numeric(y == fam_mean(family)(x, beta)),
    "logistic_deviance" = function(x, y, beta, family) log_density(family)(x, y, beta)
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


generate_df_from_formula <- function(formula, data) {
  rhs_vars <- gsub(" ", "", strsplit(deparse(rhs(formula)), "\\+(?![^()]*\\))", perl = TRUE)[[1]])

  if (identical(rhs_vars, "1")) {
    return(data.frame(rep(1, nrow(data))))
  }

  if (identical(rhs_vars, ".")) {
    return(as.data.frame(data))
  }

  new_df <- data.frame(var = character(nrow(data)))
  for (var in rhs_vars) {
    if (grepl("poly\\(.*,\\d+(\\.\\d+)?\\)", var)) {
      a <- poly(
        data[[sub("poly\\((.*),\\d+(\\.\\d+)?\\)", "\\1", var)]],
        as.double(sub("poly\\(.*,(\\d+(\\.\\d+)?)\\)", "\\1", var))
      )

      new_df <- cbind(new_df, poly(
        data[[sub("poly\\((.*),\\d+(\\.\\d+)?\\)", "\\1", var)]],
        as.double(sub("poly\\(.*,(\\d+(\\.\\d+)?)\\)", "\\1", var))
      ))
    } else {
      new_df <- cbind(new_df, parse_var_with_operator(var, data))
    }
  }
  return(new_df[-1])
}

parse_var_with_operator <- function(var, data) {
  # Split the var by '*'
  parts <- strsplit(var, "\\*")[[1]]
  parsed_parts <- vector(mode = "list", length = length(parts))

  for (i in seq_along(parts)) {
    part <- trimws(parts[i])
    # Extract variable names from the part
    vs <- regmatches(part, gregexpr("\\b[A-Za-z0-9_]+\\b", part))[[1]]
    new_part <- part

    for (v in vs) {
      # Replace variable names with data[[v]]
      new_part <- gsub(v, paste0("data[['", v, "']]"), new_part)
    }

    # Adjust the exponent notation to remove "data[['" and "']]"
    new_part <- gsub("data\\[\\['(\\d+)'\\]\\]", "\\1", new_part)

    # Append the transformed part
    parsed_parts[[i]] <- new_part
  }

  # Combine the transformed parts with '*'
  parsed_formula <- paste(parsed_parts, collapse = "*")

  return(eval(parse(text = parsed_formula)))
}

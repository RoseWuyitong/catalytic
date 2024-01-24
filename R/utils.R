log_density <- function(regression_type) {
  switch(regression_type,
    "logistic" = function(X, Y, beta, weights = FALSE, ...) {
      theta <- matrix_multiply(X, beta)
      if (isFALSE(weights)) {
        return(-sum(theta * Y - log(1 + exp(theta))))
      } else {
        return(-sum(weights * (theta * Y - log(1 + exp(theta)))))
      }
    },
    "linear" = NA
  )
}

reg_mean <- function(regression_type) {
  switch(regression_type,
    "logistic" = function(X, beta, ...) 1 / (1 + exp(-matrix_multiply(X, beta))),
    "linear" = function(X, beta, ...) matrix_multiply(X, beta),
    stop("please enter right regression type.")
  )
}

y_value <- function(regression_type) {
  switch(regression_type,
    "logistic" = function(X, beta, ...) {
      generate_data_from_dist(
        "bernoulli",
        nrow(X),
        reg_mean("logistic")(X, beta)
      )
    },
    "linear" = function(X, beta, ...) reg_mean("linear")(X, beta)
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, distribution, ...) 2 * reg_mean(distribution)(X, beta),
    "classification_error" = function(X, beta, distribution, ...) 2 * reg_mean(distribution)(X, beta),
    "logistic_deviance" = function(X, beta, ...) matrix_multiply(X, beta)
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

is_continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
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
      stats::glm(y ~ ., data = x, family = family, weights = weights)
    },
    "pseudo_mle" = if (identical(x, 1)) {
      stats::glm(y ~ 1, family = family, weights = weights, maxit = 25)
    } else {
      stats::glm(y ~ ., data = x, family = family, weights = weights, maxit = 25)
    },
    "cauchy" = if (identical(x, 1)) {
      arm::bayesglm(y ~ 1, family = family, weights = weights)
    } else {
      arm::bayesglm(y ~ ., data = x, family = family, weights = weights)
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

generate_data_from_dist <- function(generate_dist = "bernoulli",
                                    sample_size = 10,
                                    mean_ = 0,
                                    sd = 1) {
  switch(generate_dist,
    "bernoulli" = rbinom(
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
    "logarithmic" = function(x, y, beta, reg_type) {
      # cap the value between 0.0001 and 0.9999 for better calculation process
      a_cap <- pmax(0.0001, pmin(0.9999, reg_mean(data_$reg_type)(data_$obs_x, beta_tau)))
      b_cap <- pmax(0.0001, pmin(0.9999, y))
      a * log(a_cap / b_cap) + (1 - a) * log((1 - a_cap) / (1 - b_cap))
    },
    "square_error" = function(x, y, beta, reg_type) (y - reg_mean(data_$reg_type)(data_$obs_x, beta_tau))^2,
    "classification_error" = function(x, y, beta, reg_type) as.numeric(y == reg_mean(data_$reg_type)(data_$obs_x, beta_tau)),
    "logistic_deviance" = function(x, y, beta, reg_type) {
      ifelse(identical(reg_type, "logistic"),
        log_density(reg_type)(x, y, beta),
        stop("you can only use logistic_deviance for logistic regression")
      )
    }
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, reg_type) 2 * reg_mean(reg_type)(X, beta),
    "classification_error" = function(X, beta, reg_type) 2 * reg_mean(reg_type)(X, beta),
    "logistic_deviance" = function(X, beta, reg_type) {
      ifelse(identical(reg_type, "logistic"),
        matrix_multiply(X, beta),
        stop("you can only use logistic_deviance for logistic regression")
      )
    }
  )
}

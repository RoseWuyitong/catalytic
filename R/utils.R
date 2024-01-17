log_density <- function(distribution) {
  switch(distribution,
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

dist_mean <- function(distribution) {
  switch(distribution,
    "logistic" = function(X, beta, ...) 1 / (1 + exp(-matrix_multiply(X, beta))),
    "linear" = function(X, beta, ...) matrix_multiply(X, beta)
  )
}

y_value <- function(distribution) {
  switch(distribution,
    "logistic" = function(X, beta, ...) {
      generate_data(
        "bernoulli",
        nrow(X),
        dist_mean("logistic")(X, beta)
      )
    },
    "linear" = function(X, beta, ...) dist_mean("linear")(X, beta)
  )
}

lambda_function <- function(discrepancy_function) {
  switch(discrepancy_function,
    "square_error" = function(X, beta, distribution, ...) 2 * dist_mean(distribution)(X, beta),
    "classification_error" = function(X, beta, distribution, ...) 2 * dist_mean(distribution)(X, beta),
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

  return(NA)
}

#' resampling data
#'
#' @param original_data numeric or complex or logical vectors.
#' Data used for resampling.
#' @param resample_data_size integer.
#' Number of resampled output data.
#' @param resample_method str.
#' Choose from (coordinate,
#' symmetrizing,
#' smoothing,
#' flattering,
#' simple).
#' @param smoothing_noise_variance numeric.
#' The variance used for smoothing resample method only.
#'
#' @return matrix.
#' Resampled data.
#'
#' @export
#'
#' @examples
#' ## Pass a vector to resampling, and it will transform it into matrix before further process
#' resampling(1:5, resample_method = "coordinate")
#'
#' ## Pass the same input, and it will output the memerised result
#' resampling(
#'   1:5,
#'   resample_method = "coordinate"
#' ) == resampling(
#'   1:5,
#'   resample_method = "coordinate"
#' )
#'
#' ## Pass the different input, and it will output the processed result
#' resampling(1:5,
#'   resample_data_size = 5,
#'   resample_method = "coordinate"
#' )
resampling <- memoise::memoise(
  function(original_data,
           resample_data_size = 10,
           resample_method = "simple",
           smoothing_noise_variance = 1) {
    original_data <- as.matrix(original_data)
    resampled_data <- switch(resample_method,
      "coordinate" = apply(original_data,
        2,
        sample,
        size = resample_data_size,
        replace = TRUE
      ),
      "symmetrizing" = apply(
        original_data,
        2,
        function(x) {
          sample(c(x, 2 * mean(original_data) - x),
            resample_data_size,
            replace = TRUE
          )
        }
      ),
      "smoothing" = original_data + stats::rnorm(
        length(original_data),
        0,
        sqrt(smoothing_noise_variance)
      ),
      "flattering" = apply(
        original_data,
        2,
        function(x) {
          sample(c(x, sample(c(unique(original_data)), 1)),
            resample_data_size,
            replace = TRUE
          ) # todo: change
        }
      ),
      NA
    )

    return(resampled_data)
  }
)

#' generate data from distribution
#'
#' @param from_dist str.
#' The distribution for generating data, choose from
#' (bernoulli,
#' gaussian,
#' uniform)
#' @param sample_size integer.
#' The number of output sample data
#' @param dist_mean double.
#' The mean of the distribution for generating data
#' @param dist_sd double.
#' The standard deviation of the distribution for generating data
#'
#' @return matrix. sample data
#' @export
#'
#' @examples
#' ## Generate 10 data from Gaussian distribution with mean 0 and sd 1
#' generate_data_from_dist(from_dist = "gaussian")
#'
#' ## Pass the same input, and it will output the memerised result
#' generate_data_from_dist(
#'   from_dist = "gaussian"
#' ) == generate_data_from_dist(
#'   from_dist = "gaussian"
#' )
#'
#' ## Pass the different input, and it will output the processed result
#'
#' generate_data_from_dist(
#'   from_dist = "gaussian",
#'   sample_size = 5
#' )
generate_data_from_dist <- memoise::memoise(
  function(from_dist = "uniform",
           sample_size = 10,
           dist_mean = 0,
           dist_sd = 1) {
    result <- switch(from_dist,
      "bernoulli" = stats::rbinom(
        sample_size,
        1,
        meadist_meann_
      ),
      "gaussian" = stats::rnorm(
        sample_size,
        dist_mean,
        dist_sd
      ),
      "uniform" = rep(dist_mean, sample_size)
    )

    return(as.matrix(result))
  }
)

fit <- function(covariate,
                response,
                weights = FALSE,
                method = "mle",
                family = "binomial") {
  data_ <- cbind(as.data.frame(response), as.data.frame(covariate))
  names(data_)[1] <- "outcome"

  switch(method,
    "mle" = if (isFALSE(weights)) {
      glm(outcome ~ ., data = data_, family = family)
    } else {
      glm(outcome ~ ., data = data_, family = family, weights = weights)
    },
    "pseudo_mle" = if (isFALSE(weights)) {
      glm(outcome ~ ., data = data_, family = family, maxit = 25)
    } else {
      glm(outcome ~ ., data = data_, family = family, weights = weights, maxit = 25)
    },
    "cauchy" = if (isFALSE(weights)) {
      bayesglm(outcome ~ ., data = data_, family = family)
    } else {
      bayesglm(outcome ~ ., data = data_, family = family, weights = weights)
    }
  )
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

is_continue <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}

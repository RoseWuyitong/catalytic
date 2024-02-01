source("~/catalytic/R/.warning_functions.R")

#' Initialize catalytic related data
#'
#' @param observed_covariate Observed data, x.
#' @param observed_response Observed response, y.
#' @param observed_covariate_degree Degree of the observed data, default NULL.
#' @param synthetic_data_size Number of synthetic data to be generated, default 100.
#' @param synthetic_response_generation_model Simple model for generating synthetic response, default NULL(intercept).
#' @param synthetic_response_generation_family Family for generating synthetic response, default gaussian.
#' @param regression_type Regression type, default linear. Choose from ("linear", "logistic").
#'
#' @return Object contains all the information for further process.
#' @export
#'
#' @examples
#' initialization(matrix(1:9, nrow = 3),
#' c(1, 2, 3))
#'
#' initialization(matrix(1:9, nrow=3),
#' c(1, 0, 1),
#' synthetic_response_generation_family="binomial",
#' regression_type="logistic")
initialization <- function(observed_covariate,
                           observed_response,
                           observed_covariate_degree = NULL,
                           synthetic_data_size = 100,
                           synthetic_response_generation_model = NULL,
                           synthetic_response_generation_family = "gaussian",
                           regression_type = "linear") {
  data_ <- rlang::env(
    obs_x = as.data.frame(observed_covariate),
    obs_y = as.double(observed_response),
    syn_n = as.integer(synthetic_data_size),
    syn_x = NA,
    syn_y = NA,
    recd_x = data.frame(
      Type = character(),
      Process = character(),
      stringsAsFactors = FALSE
    ),
    reg_type = regression_type,
    fit_fam = synthetic_response_generation_family
  )

  .not_same_stop(
    nrow(data_$obs_x),
    length(data_$obs_y),
    "Row number of observed_covariate and observed_response not same."
  )

  data_$obs_x <- adjust_data_type(data_$obs_x)
  data_$obs_n <- nrow(data_$obs_x)
  data_$x_col_n <- ncol(data_$obs_x)

  if (is.null(observed_covariate_degree)) {
    observed_covariate_degree <- rep(1, data_$x_col_n)
  } else if (is.double(observed_covariate_degree)) {
    .not_same_stop(
      data_$x_col_n,
      length(observed_covariate_degree),
      "length of degree of matrix is not same as the column number of matrix."
    )
  } else {
    stop("Please enter degree in right format, e.g. c(1,2,1,0,...)")
  }

  data_$x_deg <- observed_covariate_degree

  # Generate synthetic covariate
  message("generating synthetic covariate")
  data_ <- do.call("initial.auto_resampling", list(data_))

  # Generate synthetic response
  message("generating synthetic response")
  x <- if (is.null(synthetic_response_generation_model)) {
    1
  } else {
    data_$obs_x[synthetic_response_generation_model]
  }
  data_$syn_y_gen_model <- fit(x, data_$obs_y, family = data_$fit_fam)
  data_$syn_y <- y_value(data_$reg_type)(data_$syn_x[synthetic_response_generation_model],
    data_$syn_y_gen_model$coefficients)

  # Combined synthetic data and observed data
  data_$comb_x <- rbind(data_$obs_x, data_$syn_x)
  data_$comb_y <- c(data_$obs_y, data_$syn_y)

  return(data_)
}

initial.auto_resampling <- function(data_) {
  data_$syn_x <- as.data.frame(apply(data_$obs_x,
    2,
    sample,
    size = data_$syn_n,
    replace = TRUE
  ))

  for (n in 1:ncol(data_$syn_x)) {
    message(paste("Column", n, "is processing"))
    process <- "Coordinate"
    # Check continue or category
    if (is_continuous(data_$obs_x[, n])) {
      data_type <- "Continuous"
      # Check order and diversity
      if (data_$x_deg[n] >= length(unique(data_$obs_x[, n])) - 1) {
        process <- c(process, "Smoothing")
        data_$syn_x[, n] <- do.call("initial.resampling", list(data_$syn_x[, n], data_$obs_x[, n], "smoothing"))
      }
      # Check skewness
      if (abs(e1071::skewness(data_$obs_x[, n])) >= 1) {
        process <- c(process, "IQR")
        data_$syn_x[, n] <- do.call("initial.resampling", list(data_$syn_x[, n], data_$obs_x[, n], "IQR"))
      }
    } else {
      data_type <- "Category"
      tbl <- table(data_$obs_x[, n])
      # Check data balance
      if (any(tbl < data_$obs_n / (2 * length(tbl)))) {
        process <- c(process, "Flattening")
        data_$syn_x[, n] <- do.call("initial.resampling", list(data_$syn_x[, n], data_$obs_x[, n], "flattening"))
      }
    }
    data_$recd_x <- rbind(data_$recd_x, data.frame(Type = data_type, Process = paste(process, collapse = "->")))
  }
  return(data_)
}



initial.resampling <- function(x, sample_x, method) {
  message(paste(method, "resampling"))
  switch(method,
    "symmetrizing" = mean(sample_x) + sample(c(-1, 1), length(x), replace = TRUE) * (x - mean(sample_x)),
    "smoothing" = x + stats::rnorm(length(x), 0, min(diff(sort(unique(sample_x)))) / 10),
    "flattening" = (function(x) {
      x[sample(length(x), length(x) / 2)] <- sample(unique(sample_x), length(x) / 2, replace = TRUE)
      return(x)
    })(x),
    "IQR" = (function(x) {
      x[sample(length(x), length(x) / 2)] <- rnorm(length(x) / 2, median(sample_x), -get_iqr(sample_x) / (2 * qnorm(1 / 4)))
      return(x)
    })(x),
    stop(paste(method, "Resampling Method Not available"))
  )
}

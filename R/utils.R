#' Multiply two matrix
#'
#' This function multiplies two matrices or data frames. The function first
#' converts data frames to matrices (if necessary), checks the dimensions,
#' and performs the matrix multiplication.
#'
#' @param matrix1 A numeric `matrix` or `data.frame`, which will be multiplied by `matrix2`.
#' @param matrix2 A numeric `matrix` or `data.frame`, which will be multiplied by `matrix1`.
#'
#' @return A numeric `matrix`, which is the product of `matrix1` and `matrix2`.
#' @export
#' @examples
#' # Create example matrices
#' matrix1 <- matrix(data = c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
#' matrix2 <- matrix(data = c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, ncol = 4)
#'
#' multiply_two_matrices(matrix1, matrix2)
#'
#' # should return the same result
#' multiply_two_matrices(matrix2, matrix1)
multiply_matrix <- function(matrix1,
                            matrix2) {
  matrix1 <- as.matrix(matrix1)
  matrix2 <- as.matrix(matrix2)

  if (ncol(matrix1) == nrow(matrix2)) {
    return(matrix1 %*% matrix2)
  } else if (nrow(matrix1) == ncol(matrix2)) {
    return(matrix2 %*% matrix1)
  } else if (ncol(matrix1) == ncol(matrix2)) {
    return(t(matrix1) %*% matrix2)
  } else if (nrow(matrix1) == nrow(matrix2)) {
    return(matrix1 %*% t(matrix2))
  }

  stop("Please note that these two matrices cannot be multiplied.",
    call. = TRUE
  )
}


log_density_binomial <- function(Y,
                                 mean,
                                 weights = 1) {
  return(
    sum(
      weights * (mean * Y - log(1 + exp(mean)))
    )
  )
}

log_density_gaussian <- function(residuals,
                                 weights = 1) {
  return(0.5 * weights * (1 + log(2 * pi * var(residuals))))
}

log_density_grad_binomial <- function(X,
                                      Y,
                                      mean,
                                      weights = 1) {
  return(
    c(
      multiply_matrix(
        matrix1 = t(weights * (Y - mean)),
        matrix2 = cbind(1, X)
      )
    )
  )
}

second_derivative <- function(X,
                              fitted_values,
                              prior_weights) {
  return((t(cbind(1, X))) %*%
    diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
    cbind(1, X))
}

mean_of_binomial <- function(X,
                             coefficients) {
  return(
    1 / (1 + exp(-multiply_two_matrices(
      matrix1 = X,
      matrix2 = coefficients
    )))
  )
}

mean_of_gaussian <- function(X,
                             coefficients) {
  return(multiply_matrix(X, coefficients))
}

correlation_matrix <- function(corr, dim) {
  matrix(corr,
    nrow = dim,
    ncol = dim
  ) +
    diag(dim) * (1 - corr)
}

standard_error <- function(X) sd(X) / sqrt(length(X))

discrepancy_of_logarithmic_error <- function(Y, est_Y) {
  Y * log(Y / est_Y) +
    (1 - Y) * log((1 - Y) / (1 - est_Y))
}

discrepancy_of_square_error <- function(Y, est_Y) (Y - est_Y)^2

discrepancy_of_classification_error <- function(Y, est_Y) as.numeric(Y == est_Y) / length(Y)

is.binomial <- function(x) x == "binomial"

is.gaussian <- function(x) x == "gaussian"

is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}


is.data_ <- function(x) inherits(x, "data_")
is.data_lmm <- function(x) inherits(x, "data_") & x$function_name == "data_initilization_lmm"
is.data_cox <- function(x) inherits(x, "data_") & x$function_name == "data_initilization_cox"

is.lmer_formula <- function(f) grepl("\\(.*|.*\\)", deparse(f))
is.surv_formula <- function(f) grepl("Surv", deparse(get_lhs(f)))

.large_number_warning <- function(num, threshold, suggestion) {
  if (any(num > threshold)) warning(suggestion, call. = FALSE)
}

.small_number_warning <- function(num, threshold, suggestion) {
  if (any(num < threshold)) warning(suggestion, call. = FALSE)
}

.small_number_stop <- function(num, threshold, suggestion) {
  if (any(num < threshold)) stop(suggestion, call. = FALSE)
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

.large_number_stop <- function(num, threshold, suggestion) {
  if (any(num > threshold)) stop(suggestion, call. = FALSE)
}

.not_true_stop <- function(statement, suggestion) {
  if (!isTRUE(statement)) stop(suggestion, call. = FALSE)
}

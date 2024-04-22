multiply_two_matrices <- function(matrix1,
                                  matrix2,
                                  intl = TRUE) {
  matrix1 <- as.matrix(matrix1)
  matrix2 <- as.matrix(matrix2)

  if (intl) {
    if (ncol(matrix1) + 1 == nrow(matrix2)) {
      return(matrix2[1] + matrix1 %*% matrix2[-1])
    } else if (
      ncol(matrix2) + 1 == nrow(matrix1)) {
      return(matrix1[1] + matrix2 %*% matrix1[-1])
    } else if (
      nrow(matrix1) + 1 == nrow(matrix2)) {
      return(matrix2[1] + matrix2[-1] %*% matrix1)
    } else if (
      nrow(matrix2) + 1 == nrow(matrix1)) {
      return(matrix1[1] + matrix1[-1] %*% matrix2)
    }
  }

  if (ncol(matrix1) == nrow(matrix2)) {
    return(matrix1 %*% matrix2)
  } else if (
    nrow(matrix1) == ncol(matrix2)) {
    return(matrix2 %*% matrix1)
  }

  stop(paste("Please note that the two matrices cannot be multiplied:",
             str(head(matrix1)),
             str(head(matrix2))))
}


log_density_binomial <- function(Y,
                                 mean,
                                 weights=1) {
  return(
    sum(
      weights * (mean * Y - log(1 + exp(mean)))
    )
  )
}

log_density_gaussian <- function(residuals,
                                 weights=1) {
  return( 0.5 * weights * (1+log(2*pi*var(residuals)))
  )
}

log_density_grad_binomial <- function(X,
                                      Y,
                                      mean,
                                      weights = 1) {
  return(
    c(
      multiply_two_matrices(matrix1 = t(weights * (Y - mean)),
                            matrix2 = cbind(1, X),
                            intl=FALSE)
    )
  )
}

second_derivative <- function(X,
                              fitted_values,
                              prior_weights
                              ) {
  return((t(cbind(1, X))) %*%
           diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
           cbind(1, X))
}

mean_of_binomial <- function(X,
                              coefficients) {
  return(
    1 / (1 + exp(-multiply_two_matrices(matrix1=X,
                                        matrix2=coefficients)))
  )
}

mean_of_gaussian <- function(X,
                              coefficients){
  return(
    multiply_two_matrices(matrix1=X,
                          matrix2=coefficients)
  )
}

correlation_matrix <- function(corr, dim) matrix(corr,
                                                 nrow = dim,
                                                 ncol = dim) +
  diag(dim) * (1 - corr)

standard_error <- function(X) sd(X) / sqrt(length(X))

discrepancy_of_logarithmic_error <- function(Y, est_Y) Y * log(Y / est_Y) +
  (1 - Y) * log((1 - Y) / (1 - est_Y))

discrepancy_of_square_error <- function(Y, est_Y) (Y-est_Y)^2

discrepancy_of_classification_error <- function(Y, est_Y) as.numeric(Y == est_Y)/length(Y)

is.binomial <- function(x) x == "binomial"

is.gaussian <- function(x) x == "gaussian"

is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}


is.data_ <- function(x) inherits(x, "data_initialization")


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




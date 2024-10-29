set.seed(1)

# ----- GLM -----

# ----- COX -----

# ----- LMM -----

mock_lmm_data <- function() {
  n <- 10
  return(data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n),
    Y = rnorm(n),
    Z1 = rnorm(n),
    Z2 = rnorm(n),
    Z3 = rnorm(n),
    group = sample(1:3, n, replace = TRUE)
  ))
}

mock_cat_lmm_initialization <- function() {
  return(cat_lmm_initialization(
    formula = Y ~ X1 + X2 + X3,
    data = mock_lmm_data(),
    x_cols = c("X1", "X2", "X3"),
    z_cols = c("Z1", "Z2", "Z3"),
    y_col = c("Y"),
    group_col = c("group")
  ))
}

##########
# ---- GLM Functions ----

#' Generate GLM Data for Testing
#'
#' This helper function generates random matrices, coefficients, and response
#' variables for Gaussian and Binomial GLMs. It also fits the models using the
#' generated data.
#'
#' @return A list containing the design matrix, coefficients, response variables,
#'         and fitted models for Gaussian and Binomial GLMs.
glm_data <- function() {
  X <- matrix(rnorm(15), ncol = 3)
  Y_gaussian <- rnorm(5)
  Y_binomial <- rbinom(5, 1, 0.5)

  coefs <- rnorm(4)

  model_gaussian <- glm(Y_gaussian ~ X, family = gaussian)
  model_binomial <- glm(Y_binomial ~ X, family = binomial)

  return(list(
    X = X,
    coefs = coefs,
    Y_gaussian = Y_gaussian,
    Y_binomial = Y_binomial,
    model_gaussian = model_gaussian,
    model_binomial = model_binomial
  ))
}

cox_data <- function() {
  X <- matrix(rnorm(15), ncol = 3)
  time <- rexp(5, rate = 0.1)
  status <- rbinom(5, 1, 0.5)

  coefs <- rnorm(3)

  model <- survival::coxph(survival::Surv(time, status) ~ X)

  return(list(
    X = X,
    coefs = coefs,
    time = time,
    status = status,
    model = model
  ))
}



#
#
# # ----- mock result -----
# mock_data_glm_gaussian <- function() {
#   return(as.data.frame(matrix(
#     data = rnorm(15),
#     ncol = 3,
#     dimnames = list(c(), c("y", "x1", "x2"))
#   )))
# }
#
# mock_data_glm_binomial <- function() {
#   data <- mock_data_glm_gaussian()
#   data[, "y"] <- as.integer(runif(data[, "y"]) > 0.5)
#   return(as.data.frame(data))
# }
#
# mock_data_initialization_glm_gaussian <- function() {
#   data_ <- data_initialization_glm(
#     formula = y ~ 1,
#     family = gaussian(),
#     data = mock_data_glm_gaussian()
#   )
#   return(data_)
# }
#
# mock_data_initialization_glm_binomial <- function() {
#   data_ <- data_initialization_glm(
#     formula = y ~ 1,
#     family = binomial(),
#     data = mock_data_glm_binomial()
#   )
#   return(data_)
# }
#
#
# # TODO: amend
# mock_data_initialization_cox <- function() {
#   data <- matrix(data = rnorm(15), ncol = 3, dimnames = list(c(), c("y", "x1", "x2")))
#   data[, "y"] <- as.integer(runif(data[, "y"]) > 0.5)
#   data_ <- data_initialization_cox(
#     formula = y ~ 1,
#     family = binomial(),
#     data = data
#   )
#
#   return(data_)
# }


# ----- Expect function -----

expect_df_equal_values <- function(actual, expected) {
  # Compare data frames ignoring attributes
  result <- all.equal(actual, expected, check.attributes = FALSE)

  # If result is TRUE, the data frames are equal
  if (is.logical(result) && result) {
    expect_true(TRUE)
  } else {
    # Otherwise, fail the test and provide details
    expect_true(FALSE, info = paste("Data frames are not equal:\n", result))
  }
}

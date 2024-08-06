library(catalytic)
source("~/catalytic/tests/testthat/helper.R")

set.seed(1)

context("misc")


# ---- GLM Functions ----

test_that("get_glm_lambda function tests", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  coefs <- c(0.5, 0.2)

  # Test for square_error
  result <- get_glm_lambda("square_error", X, coefs)
  expect_equal(result, 2 * get_glm_mean("gaussian", X, coefs))

  # Test for classification_error
  result <- get_glm_lambda("classification_error", X, coefs)
  expect_equal(result, 2 * get_glm_mean("binomial", X, coefs))

  # Test for logistic_deviance
  result <- get_glm_lambda("logistic_deviance", X, coefs)
  expect_equal(result, get_linear_predictor(X, coefs))

  # Test for invalid discrepancy_method
  expect_error(get_glm_lambda("invalid_method", X, coefs))
})

# TODO: left here
test_that("get_glm_diag_approx_cov function tests", {
  data(mtcars)
  model <- glm(mpg ~ wt + hp, data = mtcars, family = gaussian())
  X <- as.matrix(mtcars[, c("wt", "hp")])

  # Test for correct calculation
  result <- get_glm_diag_approx_cov(X, model)

  # Calculating the expected result manually
  fitted_values <- model$fitted.values
  prior_weights <- model$prior.weights
  second_derivative <- t(cbind(1, X)) %*%
    diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
    cbind(1, X)
  expected_result <- diag(chol2inv(chol(second_derivative)))

  expect_equal(result, expected_result)

  # Test with a minimal dataset
  X_minimal <- matrix(c(1, 2), nrow = 1)
  model_minimal <- glm(c(1) ~ c(1), family = binomial())
  result_minimal <- get_glm_diag_approx_cov(X_minimal, model_minimal)

  expect_true(is.numeric(result_minimal))
  expect_length(result_minimal, 1)

  # Test for invalid input: X not a matrix
  expect_error(
    get_glm_diag_approx_cov("not a matrix", model),
    "X must be a matrix"
  )

  # Test for invalid input: model not a GLM model object
  expect_error(
    get_glm_diag_approx_cov(X, "not a model"),
    "model must be a GLM model object"
  )
})


test_that("check_size_and_risk_estimate_method works correctly", {
  # Test case 1: Large dataset with parametric bootstrap method
  expect_warning(
    check_size_and_risk_estimate_method(500, "parametric_bootstrap"),
    "Due to the large size of the dataset"
  )

  # Test case 2: Small dataset with cross-validation method
  expect_warning(
    check_size_and_risk_estimate_method(25, "cross_validation"),
    "The size of the observational dataset is insufficient"
  )

  # Test case 3: Adequate size with steinian estimate method
  expect_warning(
    check_size_and_risk_estimate_method(150, "steinian_estimate"),
    NA
  )

  # Test case 4: Very large dataset with steinian estimate method
  mock_menu(2) # Simulate user input 2: "NO"

  expect_error(
    check_size_and_risk_estimate_method(500, "steinian_estimate"),
    "Process terminated by user."
  )

  # Test case 5: Very large dataset with steinian estimate method
  mock_menu(1) # Simulate user input 1: "YES"

  expect_warning(
    check_size_and_risk_estimate_method(500, "steinian_estimate"),
    "Due to the large size of the dataset"
  )

  # Test case 6: Small dataset with parametric bootstrap method
  expect_warning(
    check_size_and_risk_estimate_method(25, "parametric_bootstrap"),
    NA
  )
})

test_that("get_adjusted_data function works correctly", {
  # Test data_initialization_glm(family = gaussian)
  data_gaussian <- mock_data_initialization_glm_gaussian()
  adj_data_gaussian <- get_adjusted_data(data_gaussian, ~x1)
  expect_true(all(c("adj_x", "adj_obs_x", "adj_syn_x") %in% names(adj_data_gaussian)))
  expect_true(colnames(adj_data_gaussian$adj_x) == "x1")
  expect_df_equal_values(rbind(adj_data_gaussian$adj_obs_x, adj_data_gaussian$adj_syn_x), adj_data_gaussian$adj_x)
  expect_df_equal_values(adj_data_gaussian$adj_obs_x, as.data.frame(data_gaussian$obs_x[, "x1"]))
  expect_df_equal_values(adj_data_gaussian$adj_syn_x, as.data.frame(data_gaussian$syn_x[, "x1"]))
  expect_df_equal_values(adj_data_gaussian$adj_x, as.data.frame(data_gaussian$x[, "x1"]))

  # Test data_initialization_glm(family = binomial)
  data_binomial <- mock_data_initialization_glm_binomial()
  adj_data_binomial <- get_adjusted_data(data_binomial, ~x2)
  expect_true(all(c("adj_x", "adj_obs_x", "adj_syn_x") %in% names(adj_data_binomial)))
  expect_true(colnames(adj_data_binomial$adj_x) == "x2")
  expect_df_equal_values(rbind(adj_data_binomial$adj_obs_x, adj_data_binomial$adj_syn_x), adj_data_binomial$adj_x)
  expect_df_equal_values(adj_data_binomial$adj_obs_x, as.data.frame(data_binomial$obs_x[, "x2"]))
  expect_df_equal_values(adj_data_binomial$adj_syn_x, as.data.frame(data_binomial$syn_x[, "x2"]))
  expect_df_equal_values(adj_data_binomial$adj_x, as.data.frame(data_binomial$x[, "x2"]))

  # test data_initialization_cox
  data_cox <- mock_data_initialization_cox()
  adj_data_cox <- get_adjusted_data(data_cox, ~.)
  expect_true(all(c("adj_x", "adj_obs_x", "adj_syn_x") %in% names(adj_data_cox)))
  expect_true(colnames(adj_data_cox$adj_x) == c("x1", "x2"))
  expect_df_equal_values(rbind(adj_data_cox$adj_obs_x, adj_data_cox$adj_syn_x), adj_data_cox$adj_x)
  expect_df_equal_values(adj_data_cox$adj_obs_x, as.data.frame(data_cox$obs_x))
  expect_df_equal_values(adj_data_cox$adj_syn_x, as.data.frame(data_cox$syn_x))
  expect_df_equal_values(adj_data_cox$adj_x, as.data.frame(data_cox$x))
})

test_that("get_resampled_data function works correctly", {
  data <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = factor(c("A", "B", "A", "B", "A")),
    x3 = c(1.1, 2.2, 3.3, 4.4, 5.5)
  )

  # Test resampling functionality for numeric and factor columns
  resampled <- get_resampled_data(data, resample_size = 10)
  expect_equal(nrow(resampled$resampled_data), 10)
  expect_true(all(sapply(resampled$resampled_data, is.numeric) | sapply(resampled$resampled_data, is.factor)))

  # Test resample_only = TRUE
  resampled_only <- get_resampled_data(data, resample_size = 10, resample_only = TRUE)
  expect_null(resampled_only$resampled_data_record)
  expect_equal(nrow(resampled_only$resampled_data), 10)

  # Test if resampling record has correct structure and types
  expect_true(all(resampled$resampled_data_record$Type %in% c("Continuous", "Category")))
  expect_true(all(grepl("Coordinate", resampled$resampled_data_record$Process)))

  # Check that data degrees are considered in resampling process
  data_degree <- c(1, 2, 3)
  resampled_with_degree <- get_resampled_data(data, resample_size = 10, data_degree = data_degree)
  expect_equal(ncol(resampled_with_degree$resampled_data), 3)

  # Ensure the resampled_data is updated in the function
  expect_false(identical(resampled$resampled_data, resampled_only$resampled_data))

  # Test smoothing method applied
  not_smooth_data <- data.frame(
    x1 = rnorm(10),
    x2 = c(rep(1, 4), rep(2, 4), rep(3, 2))
  )
  x_degree <- c(1, 2)
  resampled <- get_resampled_data(not_smooth_data,
    resample_size = 20,
    data_degree = x_degree
  )

  expect_true(all(sapply(resampled$resampled_data, is.numeric)))
  expect_true(all(grepl("Smoothing", resampled$resampled_data_record$Process[2])))
  expect_true(x_degree[2] < (length(unique(resampled$resampled_data[, 2])) - 1))

  # Test deskewing method applied
  skewed_data <- data.frame(
    x1 = rnorm(10),
    x2 = c(rep(0, 8), 500, 501)
  )
  resampled <- get_resampled_data(skewed_data, resample_size = 20)
  expect_true(all(sapply(resampled$resampled_data, is.numeric)))
  expect_true(grepl("Deskewing", resampled$resampled_data_record$Process[2]))
  expect_true(abs(get_skewness(resampled$resampled_data[, 2])) < abs(get_skewness(skewed_data[, 2])))


  # Test flattening method applied
  not_flat_data <- data.frame(
    x1 = rnorm(10),
    x2 = c(rep(1, 9), 2)
  )
  resampled <- get_resampled_data(not_flat_data,
    resample_size = 20
  )

  expect_true(all(sapply(resampled$resampled_data, is.numeric)))
  expect_true(grepl("Flattening", resampled$resampled_data_record$Process[2]))
  tabl_1 <- table(not_flat_data[, 2]) / nrow(not_flat_data)
  tabl_2 <- table(resampled$resampled_data[, 2]) / nrow(resampled$resampled_data)
  expect_true(tabl_1[1] > tabl_2[1] & tabl_1[2] < tabl_2[2])
})

test_that("get_resampling_methods function works correctly", {
  x <- rnorm(10)
  original_x <- rnorm(10)

  # Test symmetrizing method
  resampled_x <- get_resampling_methods(x, original_x, method = "symmetrizing")
  expect_equal(length(resampled_x), length(x))
  expect_false(identical(resampled_x, x))
  expect_false(any(is.na(resampled_x)))

  # Test smoothing method
  resampled_x <- get_resampling_methods(x, original_x, method = "smoothing")
  expect_equal(length(resampled_x), length(x))
  expect_false(identical(resampled_x, x))
  expect_false(any(is.na(resampled_x)))

  # Test flattening method
  resampled_x <- get_resampling_methods(x, original_x, method = "flattening")
  expect_equal(length(resampled_x), length(x))
  expect_false(identical(resampled_x, x))
  expect_false(any(is.na(resampled_x)))

  # Test deskewing method
  resampled_x <- get_resampling_methods(x, original_x, method = "deskewing")
  expect_equal(length(resampled_x), length(x))
  expect_false(identical(resampled_x, x))
  expect_false(any(is.na(resampled_x)))
  expect_true(all(resampled_x >= min(original_x)))
  expect_true(all(resampled_x <= max(original_x)))
})

test_that("get_standardized_data function works correctly", {
  data <- data.frame(
    x1 = c(1, 2, NA, 4),
    x2 = c("A", "B", NA, "A")
  )

  # Test for numeric and factor conversion by using defaut replace method
  standardized_data <- get_standardized_data(data)
  expect_true(all(sapply(standardized_data, is.numeric) | sapply(standardized_data, is.factor)))
  expect_true(nrow(standardized_data) == 3)
  expect_false(anyNA(standardized_data))

  # Test for numeric and factor conversion by using mean
  standardized_data_mean <- get_standardized_data(data, mean)
  expect_true(all(sapply(standardized_data_mean, is.numeric) | sapply(standardized_data_mean, is.factor)))
  expect_true(nrow(standardized_data_mean) == 4)
  expect_false(anyNA(standardized_data_mean))
})

test_that("get_sample_data generates correct data types", {
  # Test Gaussian distribution
  data <- get_sample_data()
  expect_type(data, "double")

  # Test Binomial distribution
  binomial_data <- get_sample_data(family = "binomial", n = 20, mean = 0.5)
  expect_type(binomial_data, "integer")
  expect_true(all(binomial_data %in% c(0, 1)))
  expect_true(length(binomial_data) == 20)
})

test_that("get_family_glm function works correctly", {
  # Test with family functions and their names
  expect_equal(get_family_glm(gaussian), "gaussian")
  expect_equal(get_family_glm("gaussian"), "gaussian")
  expect_equal(get_family_glm(gaussian()), "gaussian")
  expect_equal(get_family_glm("gaussian", with_link = TRUE), "gaussian [identity]")
  expect_equal(get_family_glm(gaussian(), with_link = TRUE), "gaussian [identity]")

  # Test with other families
  expect_equal(get_family_glm(binomial), "binomial")
  expect_equal(get_family_glm("binomial"), "binomial")
  expect_equal(get_family_glm(binomial()), "binomial")
  expect_equal(get_family_glm("binomial", with_link = TRUE), "binomial [logit]")
  expect_equal(get_family_glm(binomial(), with_link = TRUE), "binomial [logit]")

  # Test with invalid family inputs
  expect_error(get_family_glm("invalid_family"), "'family' must be a family.")
  expect_error(get_family_glm(123), "'family' must be a family.")
  expect_error(get_family_glm(list()), "'family' must be a family.")
})

test_that("get_mean_glm function works correctly", {
  # Test with Gaussian family
  X <- matrix(1:6, ncol = 2)
  coefs <- c(0.5, 0.5)
  result_gaussian <- get_mean_glm(gaussian(), X, coefs)
  expected_gaussian <- get_linear_predictor(as.data.frame(X), as.vector(coefs))
  expect_equal(result_gaussian, expected_gaussian)

  # Test with Gaussian family and intercept
  coefs_with_intercept <- c(1, 0.5, 0.5)
  result_gaussian_intercept <- get_mean_glm(gaussian(), X, coefs_with_intercept)
  expected_gaussian_intercept <- get_linear_predictor(as.data.frame(X), as.vector(coefs_with_intercept))
  expect_equal(result_gaussian_intercept, expected_gaussian_intercept)

  # Test with Binomial family
  coefs_binomial <- c(0.5, 0.5)
  result_binomial <- get_mean_glm(binomial(), X, coefs_binomial)
  expected_binomial <- 1 / (1 + exp(-get_linear_predictor(as.data.frame(X), as.vector(coefs_binomial))))
  expect_equal(result_binomial, expected_binomial)

  # Test with Binomial family and intercept
  coefs_binomial_intercept <- c(1, 0.5, 0.5)
  result_binomial_intercept <- get_mean_glm(binomial(), X, coefs_binomial_intercept)
  expected_binomial_intercept <- 1 / (1 + exp(-get_linear_predictor(as.data.frame(X), as.vector(coefs_binomial_intercept))))
  expect_equal(result_binomial_intercept, expected_binomial_intercept)

  # Test with invalid family
  expect_error(get_mean_glm("invalid_family", X, coefs), "'family' must be a family.")
})

test_that("get_log_density_glm function works correctly", {
  # Test with Gaussian family
  X <- matrix(1:6, ncol = 2)
  Y <- c(1, 2, 3)
  coefs <- c(0.5, 0.5)
  result_gaussian <- get_log_density_glm(gaussian(), X, Y, coefs)
  expected_gaussian <- 0.5 * sum(1 * (1 + log(2 * pi * var(Y - get_linear_predictor(X, coefs)))))
  expect_equal(result_gaussian, expected_gaussian)

  # Test with Gaussian family and intercept
  coefs_with_intercept <- c(1, 0.5, 0.5)
  result_gaussian_intercept <- get_log_density_glm(gaussian(), X, Y, coefs_with_intercept)
  expected_gaussian_intercept <- 0.5 * sum(1 * (1 + log(2 * pi * var(Y - get_linear_predictor(X, coefs_with_intercept)))))
  expect_equal(result_gaussian_intercept, expected_gaussian_intercept)

  # Test with Binomial family
  Y_binomial <- c(0, 1, 0)
  coefs_binomial <- c(0.5, 0.5)
  result_binomial <- get_log_density_glm(binomial(), X, Y_binomial, coefs_binomial)
  expected_binomial <- sum(1 * (get_linear_predictor(X, coefs_binomial) * Y_binomial - log(1 + exp(get_linear_predictor(X, coefs_binomial)))))
  expect_equal(result_binomial, expected_binomial)

  # Test with Binomial family and intercept
  coefs_binomial_intercept <- c(1, 0.5, 0.5)
  result_binomial_intercept <- get_log_density_glm(binomial(), X, Y_binomial, coefs_binomial_intercept)
  expected_binomial_intercept <- sum(1 * (get_linear_predictor(X, coefs_binomial_intercept) * Y_binomial - log(1 + exp(get_linear_predictor(X, coefs_binomial_intercept)))))
  expect_equal(result_binomial_intercept, expected_binomial_intercept)

  # Test with invalid family
  expect_error(get_log_density_glm("invalid_family", X, Y, coefs), "'family' must be a family.")
  expect_error(get_log_density_glm(poisson(), X, Y, coefs), "Invalid GLM family.")
})

test_that("get_log_density_grad_glm function works correctly", {
  # Test with Gaussian family
  X <- matrix(1:6, ncol = 2)
  Y <- c(1, 2, 3)
  coefs <- c(0.5, 0.5)
  result_gaussian <- get_log_density_grad_glm(gaussian(), X, Y, coefs)
  expected_gaussian <- as.vector(1 * (Y - get_mean_glm(gaussian(), X, coefs)) %*% cbind(1, X))
  expect_equal(result_gaussian, expected_gaussian)

  # Test with Gaussian family and intercept
  coefs_with_intercept <- c(1, 0.5, 0.5)
  result_gaussian_intercept <- get_log_density_grad_glm(gaussian(), X, Y, coefs_with_intercept)
  expected_gaussian_intercept <- as.vector(1 * (Y - get_mean_glm(gaussian(), X, coefs_with_intercept)) %*% cbind(1, X))
  expect_equal(result_gaussian_intercept, expected_gaussian_intercept)

  # Test with Binomial family
  Y_binomial <- c(0, 1, 0)
  coefs_binomial <- c(0.5, 0.5)
  result_binomial <- get_log_density_grad_glm(binomial(), X, Y_binomial, coefs_binomial)
  expected_binomial <- as.vector(1 * (Y_binomial - get_mean_glm(binomial(), X, coefs_binomial)) %*% cbind(1, X))
  expect_equal(result_binomial, expected_binomial)

  # Test with Binomial family and intercept
  coefs_binomial_intercept <- c(1, 0.5, 0.5)
  result_binomial_intercept <- get_log_density_grad_glm(binomial(), X, Y_binomial, coefs_binomial_intercept)
  expected_binomial_intercept <- as.vector(1 * (Y_binomial - get_mean_glm(binomial(), X, coefs_binomial_intercept)) %*% cbind(1, X))
  expect_equal(result_binomial_intercept, expected_binomial_intercept)

  # Test with invalid family
  expect_error(get_log_density_grad_glm("invalid_family", X, Y, coefs), "'family' must be a family.")
  expect_error(get_log_density_grad_glm(poisson(), X, Y, coefs), "Invalid GLM family.")
})

test_that("get_formula_string function works correctly", {
  # Test with a simple formula
  f <- y ~ x1 + x2
  result_f <- get_formula_string(f)
  expect_equal(result_f, "y ~ x1 + x2")

  # Test with another simple formula
  g <- z ~ w
  result_g <- get_formula_string(g)
  expect_equal(result_g, "z ~ w")

  # Test with a formula containing spaces
  h <- a ~ b + c * d
  result_h <- get_formula_string(h)
  expect_equal(result_h, "a ~ b + c * d")
})

test_that("get_formula_lhs function works correctly", {
  # Test with a simple formula
  f <- y ~ x1 + x2
  result_f <- get_formula_lhs(f)
  expect_equal(result_f, "y")

  # Test with another simple formula
  g <- z ~ w
  result_g <- get_formula_lhs(g)
  expect_equal(result_g, "z")

  # Test with a more complex formula
  h <- a ~ b + c * d
  result_h <- get_formula_lhs(h)
  expect_equal(result_h, "a")
})

test_that("get_formula_rhs function works correctly", {
  # Test without tilde
  f <- y ~ x1 + x2
  result_f <- get_formula_rhs(f)
  expect_equal(result_f, "x1 + x2")

  # Test with tilde
  g <- z ~ w
  result_g <- get_formula_rhs(g, with_tilde = TRUE)
  expect_equal(result_g, "~ w")

  # Test with a more complex formula
  h <- a ~ b * c + d
  result_h <- get_formula_rhs(h)
  expect_equal(result_h, "b * c + d")
})

# Start writing tests
test_that("get_discrepancy computes correct 'logarithmic_error'", {
  n <- 100
  X <- model.matrix(~ predictor1 + predictor2,
    data = data.frame(predictor1 = rnorm(n), predictor2 = rnorm(n))
  )
  Y <- rbinom(n, 1, 0.5)
  coefs <- c(-1, 0.5, 0.3)

  discrepancy <- get_discrepancy("logarithmic_error", "binomial", X, Y, coefs)
  expect_length(discrepancy, 1)
  expect_true(is.numeric(discrepancy))
  expect_gt(discrepancy, 0)
})

test_that("get_discrepancy computes correct error for all methods", {
  # Test 'logarithmic_error'
  discrepancy_log_error <- get_discrepancy("logarithmic_error", "binomial", X, Y, coefs)
  expect_length(discrepancy_log_error, 1)
  expect_true(is.numeric(discrepancy_log_error))
  expect_gt(discrepancy_log_error, 0)

  # Test 'square_error'
  discrepancy_sq_error <- get_discrepancy("square_error", "binomial", X, Y, coefs)
  expect_length(discrepancy_sq_error, 1)
  expect_true(is.numeric(discrepancy_sq_error))
  expect_gt(discrepancy_sq_error, 0)

  # Test 'classification_error'
  discrepancy_clf_error <- get_discrepancy("classification_error", "binomial", X, Y, coefs)
  expect_length(discrepancy_clf_error, 1)
  expect_true(is.numeric(discrepancy_clf_error))
  expect_gt(discrepancy_clf_error, 0)
  expect_lt(discrepancy_clf_error, 1)

  # Test 'logistic_deviance'
  discrepancy_log_dev <- get_discrepancy("logistic_deviance", "binomial", X, Y, coefs)
  expect_length(discrepancy_log_dev, 1)
  expect_true(is.numeric(discrepancy_log_dev))
  expect_gt(discrepancy_log_dev, 0)
})

test_that("get_lambda computes correct lambda for all methods", {
  n <- 100
  X <- model.matrix(~ predictor1 + predictor2, data = data.frame(predictor1 = rnorm(n), predictor2 = rnorm(n)))
  coefs <- c(-1, 0.5, 0.3)

  # Test 'square_error'
  lambda_sq_error <- get_lambda("square_error", X, coefs)
  expect_length(lambda_sq_error, n)
  expect_true(is.numeric(lambda_sq_error))

  # Test 'classification_error'
  lambda_clf_error <- get_lambda("classification_error", X, coefs)
  expect_length(lambda_clf_error, n)
  expect_true(is.numeric(lambda_clf_error))

  # Test 'logistic_deviance'
  lambda_log_dev <- get_lambda("logistic_deviance", X, coefs)
  expect_length(lambda_log_dev, n)
  expect_true(is.numeric(lambda_log_dev))
})

test_that("get_iqr computes correct IQR for different intervals", {
  x <- rnorm(100)

  # Test default interval (0.5)
  iqr_default <- get_iqr(x)
  expect_true(is.numeric(iqr_default))
  expect_length(iqr_default, 1)
  expect_true(iqr_default > 0)

  # Test custom interval (0.3)
  iqr_custom <- get_iqr(x, invl = 0.3)
  expect_true(is.numeric(iqr_custom))
  expect_length(iqr_custom, 1)
  expect_true(iqr_custom > 0)

  # Test edge case where interval leads to zero IQR
  iqr_zero <- get_iqr(rep(1, 100))
  expect_true(is.numeric(iqr_zero))
  expect_length(iqr_zero, 1)
  expect_equal(iqr_zero, 0)

  # Test edge case with negative values
  x_neg <- rnorm(100, mean = -10)
  iqr_neg <- get_iqr(x_neg)
  expect_true(is.numeric(iqr_neg))
  expect_length(iqr_neg, 1)
  expect_true(iqr_neg > 0)
})


test_that("get_h_matrix calculats correctly", {
  X <- matrix(1:9, ncol = 3)
  X_star <- matrix(1:6, ncol = 3)
  tau <- 0.1

  h_matrix <- get_h_matrix(X, X_star, tau)

  # Verify dimensions of the output matrix
  expect_equal(dim(h_matrix), c(ncol(X), ncol(X)))
})

test_that("is.continuous function works correctly", {
  # Test with a numeric vector with more than two unique values
  expect_true(is.continuous(c(1, 2, 3)))

  # Test with a numeric vector with exactly two unique values
  expect_false(is.continuous(c(1, 2, 1, 2)))

  # Test with a numeric vector with only one unique value
  expect_false(is.continuous(c(1, 1, 1)))

  # Test with a non-numeric vector
  expect_false(is.continuous(c("a", "b", "c")))

  # Test with a numeric vector with no elements
  expect_false(is.continuous(numeric(0)))

  # Test with a numeric vector with NA values
  expect_false(is.continuous(c(NA, NA, NA)))

  # Test with a numeric vector with mixed NA and numbers
  expect_true(is.continuous(c(1, 2, NA, 3)))

  # Test with a non-numeric vector with more than two unique values
  expect_false(is.continuous(c("apple", "banana", "cherry")))

  # Test with an integer vector with more than two unique values
  expect_true(is.continuous(c(1L, 2L, 3L)))

  # Test with a logical vector
  expect_false(is.continuous(c(TRUE, FALSE, TRUE)))
})

test_that("is.formula function works correctly", {
  # Test with a valid formula
  expect_true(is.formula(~ x + y))
  expect_true(is.formula(y ~ x))

  # Test with a string that looks like a formula
  expect_true(is.formula("~ x + y"))
  expect_true(is.formula("y ~ x"))

  # Test with an invalid formula-like object
  expect_false(is.formula(5))
  expect_false(is.formula(NULL))
  expect_false(is.formula(NA))

  # Test with a valid formula with complex structure
  expect_true(is.formula(y ~ x + z * w))
  expect_true(is.formula(y ~ I(x^2) + log(z)))

  # Test with an expression that is not a formula
  expect_false(is.formula(expression(x + y)))

  # Test with a non-formula language object
  expect_false(is.formula(quote(x + y)))
})

test_that("is.formula_glm function correctly identifies GLM formulas", {
  # Test cases for GLM formulas
  expect_true(is.formula_glm(~ 1 + x))
  expect_true(is.formula_glm(y ~ x))
  expect_true(is.formula_glm("y ~ x"))

  # Test cases for non-GLM formulas
  expect_false(is.formula_glm(~ 1 + (1 | subject)))
  expect_false(is.formula_glm(Surv(time, status) ~ x))
})

test_that("is.formula_lmm correctly identifies LMM formulas", {
  # Test cases for LMM formulas
  expect_true(is.formula_lmm(~ 1 + (1 | subject)))

  # Test cases for non-LMM formulas
  expect_false(is.formula_lmm(y ~ x))
  expect_false(is.formula_lmm(~ 1 + x))
  expect_false(is.formula_lmm("x + y"))
})

test_that("is.formula_surv correctly identifies survival formulas", {
  # Test cases for SURV formulas
  expect_true(is.formula_surv(Surv(time, status) ~ x))

  # Test cases for non-SURV formulas
  expect_false(is.formula_surv(~ x + y))
  expect_false(is.formula_surv("x + y"))
})

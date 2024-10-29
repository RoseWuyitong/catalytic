library(catalytic)
source("~/catalytic/tests/testthat/helper.R")

set.seed(1)

context("misc")

glm_data <- glm_data
cox_data <- cox_data


# ----- GLMs -----


# ----- COX -----


# ----- LMM -----
cat_init_lmm <- mock_lmm_cat_init()
lmm_data <- mock_lmm_data()

test_that("update_lmm_variance produces expected results for edge cases", {
  # Edge case with very small variances
  result_small_variance <- update_lmm_variance(
    residual_variance = 1e-10,
    random_effect_variance = 1e-10,
    obs_z_eigenvalues = c(1, 2, 3),
    syn_z_eigenvalues = c(1, 2, 3),
    obs_adjusted_residuals = c(1, 2, 3),
    syn_adjusted_residuals = c(1, 2, 3),
    tau = 1
  )
  expect_true(is.finite(result_small_variance))

  # Check monotonicity with varying tau
  result_low_tau <- update_lmm_variance(1, 1, c(1, 1.5, 2), c(1, 1.5, 2), c(1, 2, 3), c(1, 2, 3), 1)
  result_high_tau <- update_lmm_variance(1, 1, c(1, 1.5, 2), c(1, 1.5, 2), c(1, 2, 3), c(1, 2, 3), 10)
  expect_true(result_low_tau < result_high_tau) # Higher tau should increase likelihood

  # Symmetry check with swapped eigenvalues with tau equal to length of the eigenvalues
  result_obs_syn_swapped <- update_lmm_variance(0.5, 1.2, c(1, 2, 2.5), c(1, 2, 3), c(1, 2, 3), c(1, 2, 3), 3)
  result_syn_obs_swapped <- update_lmm_variance(0.5, 1.2, c(1, 2, 3), c(1, 2, 2.5), c(1, 2, 3), c(1, 2, 3), 3)
  expect_equal(result_obs_syn_swapped, result_syn_obs_swapped)
})

test_that("validate_lmm_input checks parameters correctly", {
  # Test for valid inputs
  expect_silent(validate_lmm_input(
    cat_init = cat_init_lmm,
    tau = 1,
    residual_variance_0 = 1,
    random_effect_variance_0 = 1,
    coefs_0 = rep(1, ncol(cat_init_lmm$obs_x)),
    optimize_domain = c(0, 10),
    max_iter = 100,
    tol = 0.001,
    tau_seq = c(0, 1),
    cross_validation_fold_num = 5
  ))

  # Test for cat_init generated from cat_glm_initialization
  # Test for cat_init not generated from `cat_lmm_initialization`
  cat_init_invalid <- cat_init_lmm
  cat_init_invalid$function_name <- "wrong_initialization"
  expect_error(validate_lmm_input(
    cat_init = cat_init_invalid
  ))

  # Test for tau being negative
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tau = -1
  ))

  # Test for tau being a vector
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tau = c(1, 2)
  ))

  # Test for residual_variance_0 being zero
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    residual_variance_0 = 0
  ))

  # Test for random_effect_variance_0 being zero
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    random_effect_variance_0 = 0
  ))

  # Test for coefs_0 being too short for observation x
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    coefs_0 = c(1)
  ))

  # Test for coefs_0 being too long for observation x
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    coefs_0 = rep(1, 5)
  ))

  # Test for optimize_domain being negative
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    optimize_domain = c(-1, 0)
  ))

  # Test for optimize_domain being a single value
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    optimize_domain = 1
  ))

  # Test for optimize_domain being too long
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    optimize_domain = c(1, 2, 3)
  ))

  # Test for max_iter being too zero
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    max_iter = 0
  ))

  # Test for max_iter being a vector
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    max_iter = c(1, 2)
  ))

  # Test for tol being zero
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tol = 0
  ))

  # Test for tol being a vector
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tol = c(1, 2)
  ))


  # Test for tau_seq being negative
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tau_seq = c(-1, 1)
  ))

  # Test for tau_seq being a single value
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    tau_seq = 1
  ))

  # Test for cross_validation_fold_num being zero
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    cross_validation_fold_num = 0
  ))

  # Test for cross_validation_fold_num being a vector
  expect_error(validate_lmm_input(
    formula = Y ~ 1,
    cat_init = cat_init_lmm,
    cross_validation_fold_num = c(1, 2)
  ))
})

test_that("validate_lmm_initialization_input checks parameters correctly", {
  # Test for valid inputs
  expect_silent(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for data in matrix
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = as.matrix(lmm_data),
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for empty data
  expect_error(cat_lmm_initialization(
    formula = Y ~ X1,
    data = data.frame(),
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group"
  ))

  # Test for lmm formula
  expect_silent(validate_lmm_initialization_input(
    formula = Y ~ X1 + (1 | group),
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for surv formula
  expect_error(validate_lmm_initialization_input(
    formula = survival::Surv(time, status) ~ 1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for lack of x_cols
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for lack of z_cols
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for wrong x_cols
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3", "X4"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for wrong z_cols
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3", "Z4"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for y_col being vector
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2"),
    y_col = c("Y", "X3"),
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for y_col being wrong
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y1",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = "group",
    syn_size = NULL
  ))

  # Test for group_col being a vector
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2"),
    group_col = c("group", "Z3"),
    syn_size = NULL
  ))

  # Test for group_col being wrong
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = c("group1"),
    syn_size = NULL
  ))

  # Test for syn_size is too small
  expect_warning(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = c("group"),
    syn_size = 1
  ))

  # Test for syn_size is positive and large enough
  expect_warning(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = c("group"),
    syn_size = 50
  ))

  # Test for syn_size is zero
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = c("group"),
    syn_size = 0
  ))

  # Test for syn_size is negative
  expect_error(validate_lmm_initialization_input(
    formula = Y ~ X1,
    data = lmm_data,
    x_cols = c("X1", "X2", "X3"),
    y_col = "Y",
    z_cols = c("Z1", "Z2", "Z3"),
    group_col = c("group"),
    syn_size = -1
  ))
})

# ----- General -----


##############
# ----- General -----

test_that("is.continuous works correctly", {
  # Test with continuous numeric vector
  continuous_vector <- c(1, 2, 3, 4, 5)
  expect_true(is.continuous(continuous_vector))

  # Test with categorical vector
  categorical_vector <- c("A", "B", "A", "B")
  expect_false(is.continuous(categorical_vector))

  # Test with binary numeric vector
  binary_vector <- c(0, 1, 1, 0)
  expect_false(is.continuous(binary_vector))

  # Test with a numeric vector with only two unique values
  two_unique_vector <- c(1, 1, 2, 2)
  expect_false(is.continuous(two_unique_vector))

  # Test with numeric vector with more than two unique values
  mixed_numeric_vector <- c(1, 2, 3, 4, 5, 6)
  expect_true(is.continuous(mixed_numeric_vector))
})

test_that("get_standardized_data works correctly", {
  # Test with numeric and factor columns
  data <- data.frame(
    numeric_col = c(2, 3, NA, 4),
    factor_col = factor(c("A", "B", "A", NA)),
    stringsAsFactors = FALSE
  )

  # Test with numeric NA replacement
  standardized_data <- get_standardized_data(data, na_replace = 0)
  expect_equal(standardized_data$numeric_col, c(2, 3, 0, 4))
  expect_equal(standardized_data$factor_col, factor(c("A", "B", "A", "A")))

  # Test with mode replacement for factors
  standardized_data <- get_standardized_data(data, na_replace = mean)
  expect_equal(standardized_data$numeric_col, c(2, 2, 3, 4))

  # Test with NA omission
  data_with_na_omit <- get_standardized_data(data, na_replace = stats::na.omit)
  expect_equal(nrow(data_with_na_omit), 2)
  expect_equal(ncol(data_with_na_omit), 2)

  # Test with no NA values
  data_no_na <- data.frame(
    numeric_col = c(1, 2, 4),
    factor_col = factor(c("A", "B", "A")),
    stringsAsFactors = FALSE
  )
  standardized_data_no_na <- get_standardized_data(data_no_na)
  expect_equal(standardized_data_no_na, data_no_na)
})



# ---- GLM Functions ----

test_that("get_glm_lambda works correctly", {
  X <- glm_data$X
  coefs <- glm_data$coefs

  # Expected outputs based on the given inputs
  expected_square_error <- 2 * get_glm_mean("gaussian", X, coefs)
  expected_classification_error <- 2 * get_glm_mean("binomial", X, coefs)
  expected_logistic_deviance <- get_linear_predictor(X, coefs)

  # Test for "square_error"
  expect_equal(get_glm_lambda("square_error", X, coefs), expected_square_error)

  # Test for "classification_error"
  expect_equal(get_glm_lambda("classification_error", X, coefs), expected_classification_error)

  # Test for "logistic_deviance"
  expect_equal(get_glm_lambda("logistic_deviance", X, coefs), expected_logistic_deviance)

  # Test that an error is thrown for invalid discrepancy method
  expect_error(get_glm_lambda("invalid_method", X, coefs))
})

test_that("get_glm_diag_approx_cov works correctly", {
  X <- glm_data$X
  y <- glm_data$Y_binomial

  # Fit a GLM model
  model <- glm(y ~ X, family = binomial)

  # Compute the diagonal covariance matrix
  diag_cov <- get_glm_diag_approx_cov(X, model)

  # Manually compute the expected diagonal covariance for comparison
  hessian <- (t(cbind(1, X))) %*%
    diag(model$fitted.values * (1 - model$fitted.values) * model$prior.weights) %*%
    cbind(1, X)
  expected_diag_cov <- diag(chol2inv(chol(hessian)))

  # Test that the computed and expected values match
  expect_equal(diag_cov, expected_diag_cov)

  # Check for correct error handling with mismatched dimensions (if applicable)
  expect_error(get_glm_diag_approx_cov(matrix(1:4, ncol = 2), model))
})

test_that("get_glm_family works correctly", {
  # Test retrieving the family name
  gaussian_family_name <- get_glm_family("gaussian")
  expect_equal(gaussian_family_name, "gaussian")

  # Test retrieving the family name with the link function
  poisson_family_name_with_link <- get_glm_family("poisson", with_link = TRUE)
  expect_equal(poisson_family_name_with_link, "poisson [log]")

  # Test handling invalid family input
  expect_error(get_glm_family("invalid_family"), "'family' must be a valid GLM family.")
  expect_error(get_glm_family(function() "invalid"), "'family' must be a valid GLM family.")
})

test_that("get_glm_mean works correctly", {
  X <- glm_data$X
  coefs <- glm_data$coefs

  # Test for binomial family
  mean_binomial <- get_glm_mean("binomial", X, coefs)
  expected_binomial <- 1 / (1 + exp(-get_linear_predictor(X, coefs)))
  expect_equal(mean_binomial, expected_binomial, tolerance = 1e-6)

  # Test for gaussian family
  mean_gaussian <- get_glm_mean("gaussian", X, coefs)
  expected_gaussian <- get_linear_predictor(X, coefs)
  expect_equal(mean_gaussian, expected_gaussian, tolerance = 1e-6)

  # Test for invalid family
  expect_error(get_glm_mean("invalid_family", X, coefs), "'family' must be a valid GLM family.")
})

test_that("get_glm_log_density works correctly", {
  X <- glm_data$X
  Y_binomial <- glm_data$Y_binomial
  Y_gaussian <- glm_data$Y_gaussian
  coefs <- glm_data$coefs
  weights <- rnorm(length(Y_binomial)) # Using radom weights

  # Test for binomial family
  log_density_binomial <- get_glm_log_density("binomial", X, Y_binomial, coefs, weights)
  lp_binomial <- get_linear_predictor(X, coefs)
  expected_log_density_binomial <- sum(weights * (Y_binomial * lp_binomial - log(1 + exp(lp_binomial))))
  expect_equal(log_density_binomial, expected_log_density_binomial, tolerance = 1e-6)

  # Test for gaussian family
  log_density_gaussian <- get_glm_log_density("gaussian", X, Y_gaussian, coefs, weights)
  residuals <- Y_gaussian - get_linear_predictor(X, coefs)
  expected_log_density_gaussian <- 0.5 * sum(weights * (log(2 * pi * var(residuals)) + 1))
  expect_equal(log_density_gaussian, expected_log_density_gaussian, tolerance = 1e-6)

  # Test for invalid family
  expect_error(get_glm_log_density("invalid_family", X, Y_binomial, coefs), "'family' must be a valid GLM family.")
})

test_that("get_glm_log_density_grad works correctly", {
  X <- glm_data$X
  Y_binomial <- glm_data$Y_binomial
  Y_gaussian <- glm_data$Y_gaussian
  coefs <- glm_data$coefs

  # Retrieve GLM models from glm_data
  model_binomial <- glm_data$model_binomial
  model_gaussian <- glm_data$model_gaussian

  # Compute gradients
  grad_binomial <- get_glm_log_density_grad("binomial", X, Y_binomial, coefs)
  grad_gaussian <- get_glm_log_density_grad("gaussian", X, Y_gaussian, coefs)

  # Manually compute expected gradients for comparison
  mean_binomial <- get_glm_mean("binomial", X, coefs)
  grad_binomial_expected <- (model_binomial$prior.weights * (Y_binomial - mean_binomial)) %*% cbind(1, X)

  mean_gaussian <- get_glm_mean("gaussian", X, coefs)
  grad_gaussian_expected <- (model_gaussian$prior.weights * (Y_gaussian - mean_gaussian)) %*% cbind(1, X)

  # Test that computed values match expected values
  expect_equal(grad_binomial, grad_binomial_expected[1, ])
  expect_equal(grad_gaussian, grad_gaussian_expected[1, ])
})

test_that("get_glm_sample_data works correctly", {
  # Test data generation for binomial family
  binomial_data <- get_glm_sample_data(family = "binomial", n = 20, mean = 0.3)
  expect_equal(length(binomial_data), 20)
  expect_true(all(binomial_data %in% c(0, 1)))

  # Test data generation for Gaussian family
  gaussian_data <- get_glm_sample_data(family = "gaussian", n = 20, mean = 5, sd = 2)
  expect_equal(length(gaussian_data), 20)
  expect_true(all(is.numeric(gaussian_data)))

  # Check for error with an invalid family
  expect_error(get_glm_sample_data(family = "invalid_family", n = 20))
})

test_that("check_glm_tune_size_and_risk_estimate_method behaves correctly", {
  # Mock function to replace 'menu' for testing
  mock_menu <- function(choices, title) {
    return(1) # Simulates user choosing "YES"
  }


  # Mock the menu function
  assign("menu", mock_menu, envir = globalenv())

  # Test for large dataset with 'parametric_bootstrap' method
  expect_warning(
    check_glm_tune_size_and_risk_estimate_method(obs_size = 700, risk_estimate_method = "parametric_bootstrap"),
    "Due to the large size of the dataset, Please consider utilizing `risk_estimate_method = cross_validation` for quicker estimate process."
  )

  # Test for large dataset with 'steinian_estimate' method
  expect_warning(
    check_glm_tune_size_and_risk_estimate_method(obs_size = 500, risk_estimate_method = "steinian_estimate"),
    "Due to the large size of the dataset, Please consider utilizing `risk_estimate_method = cross_validation` for quicker estimate process."
  )

  # Test for small dataset with 'cross_validation' method
  expect_warning(
    check_glm_tune_size_and_risk_estimate_method(obs_size = 20, risk_estimate_method = "cross_validation"),
    "The size of the observational dataset is insufficient for cross validation to yield optimal results."
  )

  # Test for medium-sized dataset with 'parametric_bootstrap' method (no warning)
  expect_silent(
    check_glm_tune_size_and_risk_estimate_method(obs_size = 250, risk_estimate_method = "parametric_bootstrap")
  )
})


# ---- COX Functions ----

test_that("get_cox_qr_solve works correctly", {
  # Define valid input
  coef_matrix <- matrix(c(2, 1, 1, 2), ncol = 2)
  constant_vector <- c(5, 6)

  # Expected solution
  expected_solution <- solve(coef_matrix, constant_vector)

  # Compute the solution
  solution <- get_cox_qr_solve(coef_matrix, constant_vector)

  # Test that the computed solution matches the expected solution
  expect_equal(solution, expected_solution)

  # Test with a matrix that might produce NA values
  coef_matrix_na <- matrix(c(1, 1, 1, 1), ncol = 2)
  constant_vector_na <- c(1, 1)

  # Compute the solution
  solution_na <- get_cox_qr_solve(coef_matrix_na, constant_vector_na)

  # Test that the solution does not contain NA values and is a vector
  expect_true(is.numeric(solution_na))
  expect_false(any(is.na(solution_na)))
  expect_equal(length(solution_na), length(constant_vector_na))

  # Test with another case that produces errors
  coef_matrix_error <- matrix(c(0, 0, 0, 0), ncol = 2)
  constant_vector_error <- c(1, 1)

  # Compute the solution
  solution_error <- get_cox_qr_solve(coef_matrix_error, constant_vector_error)

  # Test that the solution is a vector of default values
  expect_true(all(solution_error == 0.0001))
  expect_equal(length(solution_error), length(constant_vector_error))
})

test_that("get_cox_hessian works correctly", {
  X <- cox_data$X
  time <- cox_data$time
  status <- cox_data$status
  coefs <- cox_data$coefs

  # Compute Hessian matrix
  hessian_matrix <- get_cox_hessian(X, time, status, coefs)

  # Test if the computed Hessian is a matrix
  expect_true(is.matrix(hessian_matrix))

  # Test the dimensions of the Hessian matrix
  expect_equal(dim(hessian_matrix), c(ncol(X), ncol(X)))

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  status_small <- c(1)
  coefs_small <- c(1)
  expect_silent(get_cox_hessian(X_small, time_small, status_small, coefs_small))

  # Test with entry points provided
  entry_points <- c(0, 0, 0, 0, 0)
  expect_silent(get_cox_hessian(X, time, status, coefs, entry_points))

  # Test with no events
  status_no_event <- rbinom(5, 1, 0)
  expect_silent(get_cox_hessian(X, time, status_no_event, coefs))
})

test_that("get_cox_syn_hessian works correctly", {
  # Define test data from cox_data
  X <- cox_data$X
  time <- cox_data$time
  coefs <- cox_data$coefs
  hazard_constant <- 0.7

  # Compute Hessian using the function
  hessian <- get_cox_syn_hessian(X, time, coefs, hazard_constant)

  # Test that the Hessian is a numeric matrix
  expect_true(is.matrix(hessian))
  expect_true(is.numeric(hessian))

  # Test that the Hessian has the correct dimensions (number of columns in X)
  expect_equal(dim(hessian), c(ncol(X), ncol(X)))

  # Test that the Hessian is symmetric (should be if correctly calculated)
  expect_true(all.equal(hessian, t(hessian), tolerance = 1e-8))

  # Test with minimal input (edge case)
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  coefs_small <- c(1)
  hazard_constant_small <- 0.5

  # Compute the Hessian with minimal input
  hessian_small <- get_cox_syn_hessian(X_small, time_small, coefs_small, hazard_constant_small)

  # Test that the Hessian is a 1x1 numeric matrix
  expect_true(is.matrix(hessian_small))
  expect_true(is.numeric(hessian_small))
  expect_equal(dim(hessian_small), c(1, 1))

  # Test that the Hessian with minimal input has a valid value (should be a negative number)
  expect_lt(hessian_small, 0)

  # Test with no variation in time (all times the same)
  time_same <- rep(1, nrow(X))
  hessian_time_same <- get_cox_syn_hessian(X, time_same, coefs, hazard_constant)

  # Test that the Hessian is still a valid numeric matrix
  expect_true(is.matrix(hessian_time_same))
  expect_true(is.numeric(hessian_time_same))
  expect_equal(dim(hessian_time_same), c(ncol(X), ncol(X)))

  # Test with a different hazard constant
  hazard_constant_diff <- 1.2
  hessian_diff <- get_cox_syn_hessian(X, time, coefs, hazard_constant_diff)

  # Test that the Hessian changes with a different hazard constant
  expect_false(isTRUE(hessian == hessian_diff))
})


test_that("get_cox_gradient works correctly", {
  # Define test data
  X <- cox_data$X
  time <- cox_data$time
  status <- cox_data$status
  coefs <- cox_data$coefs

  # Compute gradient
  gradient <- get_cox_gradient(X, time, status, coefs)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient))

  # Test that the length of the gradient matches the number of columns in X
  expect_equal(length(gradient), ncol(X))

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  status_small <- c(1)
  coefs_small <- c(1)

  # Compute the gradient with minimal input
  gradient_small <- get_cox_gradient(X_small, time_small, status_small, coefs_small)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_small))
  expect_equal(length(gradient_small), ncol(X_small))

  # Test with entry points provided
  entry_points <- c(0, 0, 0, 0, 0)
  gradient_entry <- get_cox_gradient(X, time, status, coefs, entry_points)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_entry))
  expect_equal(length(gradient_entry), ncol(X))

  # Test with no events
  status_no_event <- rep(0, length(status))
  gradient_no_event <- get_cox_gradient(X, time, status_no_event, coefs)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_no_event))
  expect_equal(length(gradient_no_event), ncol(X))
})

test_that("get_cox_syn_gradient works correctly", {
  X <- cox_data$X
  time <- cox_data$time
  coefs <- cox_data$coefs
  hazard_constant <- 0.5 # Example hazard constant

  # Compute gradient
  gradient <- get_cox_syn_gradient(X, time, coefs, hazard_constant)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient))

  # Test that the length of the gradient matches the number of columns in X
  expect_equal(length(gradient), ncol(X))

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  coefs_small <- c(1)
  hazard_constant_small <- 0.5

  # Compute the gradient with minimal input
  gradient_small <- get_cox_syn_gradient(X_small, time_small, coefs_small, hazard_constant_small)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_small))
  expect_equal(length(gradient_small), ncol(X_small))

  # Test with a different hazard constant
  hazard_constant_diff <- 1.0
  gradient_diff <- get_cox_syn_gradient(X, time, coefs, hazard_constant_diff)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_diff))
  expect_equal(length(gradient_diff), ncol(X))
})
test_that("get_cox_syn_gradient works correctly", {
  X <- cox_data$X
  time <- cox_data$time
  coefs <- cox_data$coefs
  hazard_constant <- 0.5 # Example hazard constant

  # Compute gradient
  gradient <- get_cox_syn_gradient(X, time, coefs, hazard_constant)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient))

  # Test that the length of the gradient matches the number of columns in X
  expect_equal(length(gradient), ncol(X))

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  coefs_small <- c(1)
  hazard_constant_small <- 0.5

  # Compute the gradient with minimal input
  gradient_small <- get_cox_syn_gradient(X_small, time_small, coefs_small, hazard_constant_small)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_small))
  expect_equal(length(gradient_small), ncol(X_small))

  # Test with a different hazard constant
  hazard_constant_diff <- 1.0
  gradient_diff <- get_cox_syn_gradient(X, time, coefs, hazard_constant_diff)

  # Test that the gradient is a numeric vector
  expect_true(is.numeric(gradient_diff))
  expect_equal(length(gradient_diff), ncol(X))
})

test_that("get_cox_partial_likelihood works correctly", {
  X <- cox_data$X
  time <- cox_data$time
  status <- cox_data$status
  coefs <- cox_data$coefs
  entry_points <- rep(0, length(time)) # Optional, can be omitted

  # Compute partial likelihood
  partial_likelihood <- get_cox_partial_likelihood(X, time, status, coefs, entry_points)

  # Test that the partial likelihood is a numeric scalar
  expect_true(is.numeric(partial_likelihood))
  expect_length(partial_likelihood, 1)

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  time_small <- c(1)
  status_small <- c(1)
  coefs_small <- c(1)

  # Compute the partial likelihood with minimal input
  partial_likelihood_small <- get_cox_partial_likelihood(X_small, time_small, status_small, coefs_small)

  # Test that the partial likelihood is a numeric scalar
  expect_true(is.numeric(partial_likelihood_small))
  expect_length(partial_likelihood_small, 1)
})

test_that("get_cox_kappa works correctly", {
  X <- cox_data$X
  time <- cox_data$time
  status <- cox_data$status
  coefs <- cox_data$coefs

  # Estimate kappa
  kappa <- get_cox_kappa(X, time, status)

  # Test that kappa is a numeric value
  expect_true(is.numeric(kappa))

  # Test that kappa is a scalar
  expect_equal(length(kappa), 1)

  # Test with initial values
  initial_values <- rep(0.1, ncol(X))
  kappa_initial <- get_cox_kappa(X, time, status, initial_values = initial_values)
  expect_true(is.numeric(kappa_initial))
  expect_equal(length(kappa_initial), 1)

  # Test with entry points provided
  entry_points <- c(0, 0, 0, 0, 0)
  kappa_entry <- get_cox_kappa(X, time, status, entry_points = entry_points)
  expect_true(is.numeric(kappa_entry))
  expect_equal(length(kappa_entry), 1)

  # Test with a small tolerance to force early stopping
  kappa_tol <- get_cox_kappa(X, time, status, tol = 1e-10)
  expect_true(is.numeric(kappa_tol))
  expect_equal(length(kappa_tol), 1)
})

test_that("get_cox_risk_and_failure_sets works correctly", {
  time_vector <- cox_data$time
  status_vector <- cox_data$status
  indicator_vector <- sort(unique(c(time_vector, time_vector + rexp(length(time_vector)))))

  # Calculate risk and failure sets
  sets <- get_cox_risk_and_failure_sets(time_vector, status_vector, indicator_vector)

  # Test that the result is a list with two elements
  expect_true(is.list(sets))
  expect_equal(length(sets), 2)

  # Test that the risk_set and failure_set are matrices
  expect_true(is.matrix(sets$risk_set))
  expect_true(is.matrix(sets$failure_set))

  # Test that the dimensions of risk_set and failure_set match the input
  expect_equal(dim(sets$risk_set), c(length(time_vector), length(indicator_vector)))
  expect_equal(dim(sets$failure_set), c(length(time_vector), length(indicator_vector)))

  # Test that failure_set has only one 1 per row where status = 1
  expect_true(all(rowSums(sets$failure_set[status_vector == 1, ]) == 1))

  # Test that failure_set is all zeros where status = 0
  expect_true(all(sets$failure_set[status_vector == 0, ] == 0))

  # Test that risk_set is all zeros or ones
  expect_true(all(sets$risk_set %in% c(0, 1)))
})

test_that("get_cox_risk_set_idx works correctly", {
  # Define test data from cox_data
  time_of_interest <- 5
  entry_vector <- cox_data$entry
  time_vector <- cox_data$time
  status_vector <- cox_data$status

  # Compute risk set indices
  risk_set_idx <- get_cox_risk_set_idx(time_of_interest, entry_vector, time_vector, status_vector)

  # Test that the output is a vector
  expect_true(is.vector(risk_set_idx))

  # Test that all indices are within the range of the time_vector length
  expect_true(all(risk_set_idx >= 1 & risk_set_idx <= length(time_vector)))

  # Test that all subjects in the risk set entered before or at the time of interest
  expect_true(all(time_of_interest >= entry_vector[risk_set_idx]))

  # Test that all subjects in the risk set have either the event at the time of interest or are still at risk
  expect_true(all((time_vector[risk_set_idx] == time_of_interest & status_vector[risk_set_idx] == 1) |
    (time_vector[risk_set_idx] + 1e-08 > time_of_interest)))

  # Test with a time_of_interest that is not in time_vector
  time_of_interest_new <- max(time_vector) + 1
  risk_set_idx_new <- get_cox_risk_set_idx(time_of_interest_new, entry_vector, time_vector, status_vector)
  expect_true(length(risk_set_idx_new) == 0)
})

# ---- LMM Functions ----

test_that("update_lmm_sigma_square_cat works correctly", {
  sigma_square <- 1
  sigma_u_square <- 0.5
  d_seq <- c(2, 3)
  d_seq_star <- c(1, 2)
  q <- c(0.5, 1)
  q_star <- c(0.3, 0.7)
  tau <- 0.1

  # Compute log marginal likelihood
  log_marginal_likelihood <- update_lmm_sigma_square_cat(
    sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau
  )

  # Test that the output is a numeric scalar
  expect_true(is.numeric(log_marginal_likelihood))
  expect_length(log_marginal_likelihood, 1)

  # Test with minimal input
  sigma_square_small <- 0.1
  sigma_u_square_small <- 0.2
  d_seq_small <- c(1)
  d_seq_star_small <- c(1)
  q_small <- c(0.1)
  q_star_small <- c(0.1)
  tau_small <- 0.01

  log_marginal_likelihood_small <- update_lmm_sigma_square_cat(
    sigma_square_small, sigma_u_square_small, d_seq_small, d_seq_star_small, q_small, q_star_small, tau_small
  )

  expect_true(is.numeric(log_marginal_likelihood_small))
  expect_length(log_marginal_likelihood_small, 1)

  # Test with a larger tau
  tau_large <- 1.0
  log_marginal_likelihood_large_tau <- update_lmm_sigma_square_cat(
    sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau_large
  )

  expect_true(is.numeric(log_marginal_likelihood_large_tau))
  expect_length(log_marginal_likelihood_large_tau, 1)

  expect_true(is.numeric(log_marginal_likelihood_zero_sigma))
  expect_length(log_marginal_likelihood_zero_sigma, 1)
})

# ---- Other Functions ----

test_that("get_h_matrix works correctly", {
  # Define test data
  obs_x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  syn_x <- matrix(c(7, 8, 9, 10, 11, 12), ncol = 2)
  tau <- 0.1

  # Compute H matrix
  H_matrix <- get_h_matrix(obs_x, syn_x, tau)

  # Check that the output is a matrix
  expect_true(is.matrix(H_matrix))

  # Test dimensions of the result
  expected_nrow <- nrow(obs_x)
  expected_ncol <- nrow(obs_x)
  expect_equal(dim(H_matrix), c(expected_nrow, expected_ncol))

  # Test with minimal input
  obs_x_small <- matrix(c(1), ncol = 1)
  syn_x_small <- matrix(c(2), ncol = 1)
  tau_small <- 0.01

  H_matrix_small <- get_h_matrix(obs_x_small, syn_x_small, tau_small)
  expect_true(is.matrix(H_matrix_small))
  expect_equal(dim(H_matrix_small), c(1, 1))
})

test_that("get_adjusted_data works correctly", {
  # Create a sample cat_init object
  cat_init <- list(
    x = data.frame(a = c(1, 2, 3, 4, 5), b = c(6, 7, 8, 9, 10)),
    obs_size = 3,
    size = 5
  )
  formula <- "x ~ b"

  # Adjust the data
  adjusted_data <- get_adjusted_data(cat_init, formula)

  # Test if adjusted_data contains the expected components
  expect_true("adj_x" %in% names(adjusted_data))
  expect_true("adj_obs_x" %in% names(adjusted_data))
  expect_true("adj_syn_x" %in% names(adjusted_data))

  # Test the dimensions of the adjusted data
  expect_equal(nrow(adjusted_data$adj_x), 5)
  expect_equal(ncol(adjusted_data$adj_x), 1)

  expect_equal(nrow(adjusted_data$adj_obs_x), 3)
  expect_equal(nrow(adjusted_data$adj_syn_x), 2)

  # Test if the adjusted data has numeric columns
  expect_true(all(sapply(adjusted_data$adj_x, is.numeric)))
  expect_true(all(sapply(adjusted_data$adj_obs_x, is.numeric)))
  expect_true(all(sapply(adjusted_data$adj_syn_x, is.numeric)))

  # Test column names consistency
  expect_equal(colnames(adjusted_data$adj_obs_x), colnames(adjusted_data$adj_syn_x))
  expect_equal(colnames(adjusted_data$adj_obs_x), colnames(adjusted_data$adj_x))
})

test_that("get_resampling_methods works correctly", {
  # Create sample data
  original_x <- rnorm(100)
  x <- original_x + rnorm(100, 0, 0.5)

  # Test the "symmetrizing" method
  sym_x <- get_resampling_methods(x, original_x, method = "symmetrizing")
  expect_equal(length(sym_x), length(x))

  # Test the "smoothing" method
  smoothed_x <- get_resampling_methods(x, original_x, method = "smoothing")
  expect_equal(length(smoothed_x), length(x))

  # Test the "flattening" method
  flattened_x <- get_resampling_methods(x, original_x, method = "flattening")
  expect_equal(length(flattened_x), length(x))

  # Test the "deskewing" method
  deskewed_x <- get_resampling_methods(x, original_x, method = "deskewing")
  expect_equal(length(deskewed_x), length(x))

  # Check if all methods return numeric vectors
  expect_true(is.numeric(sym_x))
  expect_true(is.numeric(smoothed_x))
  expect_true(is.numeric(flattened_x))
  expect_true(is.numeric(deskewed_x))
})


test_that("get_linear_predictor works correctly", {
  # Test with no intercept
  X <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  coefs <- c(0.5, 1.5)

  linear_predictor <- get_linear_predictor(X, coefs)
  expected_predictor <- c(0.5 * 1 + 1.5 * 3, 0.5 * 2 + 1.5 * 4)
  expect_equal(linear_predictor, expected_predictor)

  # Test with intercept
  coefs_with_intercept <- c(1, 0.5, 1.5)

  linear_predictor_with_intercept <- get_linear_predictor(X, coefs_with_intercept)
  expected_predictor_with_intercept <- c(1 + 0.5 * 1 + 1.5 * 3, 1 + 0.5 * 2 + 1.5 * 4)
  expect_equal(linear_predictor_with_intercept, expected_predictor_with_intercept)

  # Test with dimension mismatch
  expect_error(
    get_linear_predictor(X, c(0.5)),
    "Dimensions of coefs and predictor variables do not match"
  )

  # Test with minimal input
  X_small <- matrix(1, ncol = 1)
  coefs_small <- c(1)
  expect_equal(get_linear_predictor(X_small, coefs_small), 1)

  # Test with intercept and minimal input
  coefs_small_with_intercept <- c(1, 2)
  expect_equal(get_linear_predictor(X_small, coefs_small_with_intercept), 3)
})

test_that("get_formula_string works correctly", {
  # Test with a simple formula
  formula <- y ~ x1 + x2
  expected_string <- "y ~ x1 + x2"
  expect_equal(get_formula_string(formula), expected_string)

  # Test with an interaction term
  formula_interaction <- y ~ x1 * x2
  expected_string_interaction <- "y ~ x1 * x2"
  expect_equal(get_formula_string(formula_interaction), expected_string_interaction)

  # Test with a more complex formula
  formula_complex <- y ~ x1 + x2 + I(x1^2) + log(x2)
  expected_string_complex <- "y ~ x1 + x2 + I(x1^2) + log(x2)"
  expect_equal(get_formula_string(formula_complex), expected_string_complex)
})

test_that("get_formula_lhs works correctly", {
  # Test with a simple formula
  formula <- y ~ x1 + x2
  expected_lhs <- "y"
  expect_equal(get_formula_lhs(formula), expected_lhs)

  # Test with an interaction term
  formula_interaction <- y ~ x1 * x2
  expected_lhs_interaction <- "y"
  expect_equal(get_formula_lhs(formula_interaction), expected_lhs_interaction)

  # Test with a more complex formula
  formula_complex <- y ~ x1 + x2 + I(x1^2) + log(x2)
  expected_lhs_complex <- "y"
  expect_equal(get_formula_lhs(formula_complex), expected_lhs_complex)
})

test_that("get_formula_rhs works correctly", {
  # Test with a simple formula
  formula <- y ~ x1 + x2
  expected_rhs <- "x1 + x2"
  expect_equal(get_formula_rhs(formula), expected_rhs)

  # Test with tilde included
  expected_rhs_with_tilde <- "~ x1 + x2"
  expect_equal(get_formula_rhs(formula, with_tilde = TRUE), expected_rhs_with_tilde)

  # Test with interaction term
  formula_interaction <- y ~ x1 * x2
  expected_rhs_interaction <- "x1 * x2"
  expect_equal(get_formula_rhs(formula_interaction), expected_rhs_interaction)

  # Test with a more complex formula
  formula_complex <- y ~ x1 + x2 + I(x1^2) + log(x2)
  expected_rhs_complex <- "x1 + x2 + I(x1^2) + log(x2)"
  expect_equal(get_formula_rhs(formula_complex), expected_rhs_complex)

  # Test with tilde for complex formula
  expected_rhs_complex_with_tilde <- "~ x1 + x2 + I(x1^2) + log(x2)"
  expect_equal(get_formula_rhs(formula_complex, with_tilde = TRUE), expected_rhs_complex_with_tilde)
})

test_that("get_discrepancy works correctly", {
  # Sample data for testing
  X <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  coefs <- c(0.5, 1.5)
  Y <- c(1, 0)
  est_Y <- get_glm_mean("binomial", X, coefs)

  # Test logarithmic error
  discrepancy_log <- get_discrepancy("logarithmic_error", X = X, Y = Y, coefs = coefs, est_Y = est_Y)
  expected_log <- mean(Y * log(pmax(0.0001, pmin(0.9999, Y)) / pmax(0.0001, pmin(0.9999, est_Y))) +
    (1 - Y) * log((1 - pmax(0.0001, pmin(0.9999, Y))) / (1 - pmax(0.0001, pmin(0.9999, est_Y)))))
  expect_equal(discrepancy_log, expected_log, tolerance = 0.01)

  # Test square error
  discrepancy_sq <- get_discrepancy("square_error", X = X, Y = Y, coefs = coefs, est_Y = est_Y)
  expected_sq <- mean((Y - est_Y)^2)
  expect_equal(discrepancy_sq, expected_sq, tolerance = 1e-6)

  # Test classification error
  discrepancy_class <- get_discrepancy("classification_error", X = X, Y = Y, coefs = coefs, est_Y = est_Y)
  expected_class <- mean(Y - 2 * Y * est_Y + est_Y)
  expect_equal(discrepancy_class, expected_class, tolerance = 1e-6)

  # Test logistic deviance
  # Note: Requires a valid implementation of get_glm_log_density
  discrepancy_logistic <- get_discrepancy("logistic_deviance", family = "binomial", X = X, Y = Y, coefs = coefs)
  expected_logistic <- -get_glm_log_density("binomial", X, Y, coefs)
  expect_equal(discrepancy_logistic, expected_logistic, tolerance = 1e-6)
})

# ---- Identification Functions ----


test_that("is.formula works correctly", {
  # Test with a formula object
  formula_obj <- y ~ x + z
  expect_true(is.formula(formula_obj))

  # Test with a non-formula object (numeric vector)
  non_formula_obj <- c(1, 2, 3)
  expect_false(is.formula(non_formula_obj))

  # Test with a non-formula object (character string)
  formula_string <- "y ~ x + z"
  expect_true(is.formula(formula_string))

  # Test with a formula that has no tilde (not really a formula)
  no_tilde_obj <- "y + x"
  expect_false(is.formula(no_tilde_obj))
})

test_that("is.formula_glm works correctly", {
  # Mock functions for formula checks
  is.formula_lmm <- function(f) {
    grepl("\\(1\\|", deparse(f))
  }
  is.formula_surv <- function(f) {
    grepl("Surv\\(", deparse(f))
  }

  # Test with a GLM formula
  glm_formula <- y ~ x + z
  expect_true(is.formula_glm(glm_formula))

  # Test with an LMM formula
  lmm_formula <- y ~ x + (1 | group)
  expect_false(is.formula_glm(lmm_formula))

  # Test with a Survival Analysis formula
  surv_formula <- Surv(time, status) ~ x
  expect_false(is.formula_glm(surv_formula))

  # Test with a non-formula object
  non_formula_obj <- c(1, 2, 3)
  expect_false(is.formula_glm(non_formula_obj))
})

test_that("is.formula_lmm works correctly", {
  # Example with an LMM formula
  formula_lmm <- as.formula("response ~ predictor + (1 | group)")
  expect_true(is.formula_lmm(formula_lmm))

  # Example with a GLM formula
  formula_glm <- as.formula("response ~ predictor")
  expect_false(is.formula_lmm(formula_glm))

  # Example with a formula for a survival model
  formula_surv <- as.formula("Surv(time, status) ~ predictor")
  expect_false(is.formula_lmm(formula_surv))

  # Example with a non-formula object
  non_formula <- "response ~ predictor"
  expect_false(is.formula_lmm(non_formula))
})

test_that("is.formula_surv works correctly", {
  # Example with a survival formula
  formula_surv <- as.formula("Surv(time, status) ~ predictor")
  expect_true(is.formula_surv(formula_surv))

  # Example with a GLM formula
  formula_glm <- as.formula("response ~ predictor")
  expect_false(is.formula_surv(formula_glm))

  # Example with an LMM formula
  formula_lmm <- as.formula("response ~ predictor + (1 | group)")
  expect_false(is.formula_surv(formula_lmm))

  # Example with a non-formula object
  non_formula <- "Surv(time, status) ~ predictor"
  expect_true(is.formula_surv(non_formula))
})

test_that("is.cat_glm_init works correctly", {
  # Example of creating a cat_initialization object with correct function_name
  init_obj <- list(function_name = "cat_glm_initialization")
  class(init_obj) <- "cat_initialization"
  expect_true(is.cat_glm_init(init_obj))

  # Example with a different function_name
  init_obj_diff <- list(function_name = "cat_other_initialization")
  class(init_obj_diff) <- "cat_initialization"
  expect_false(is.cat_glm_init(init_obj_diff))

  # Example with a non-cat_initialization object
  non_cat_obj <- list(function_name = "cat_glm_initialization")
  expect_false(is.cat_glm_init(non_cat_obj))
})

test_that("is.cat_lmm_init works correctly", {
  # Example of creating a cat_initialization object with correct function_name
  init_obj <- list(function_name = "cat_lmm_initialization")
  class(init_obj) <- "cat_initialization"
  expect_true(is.cat_lmm_init(init_obj))

  # Example with a different function_name
  init_obj_diff <- list(function_name = "cat_other_initialization")
  class(init_obj_diff) <- "cat_initialization"
  expect_false(is.cat_lmm_init(init_obj_diff))

  # Example with a non-cat_initialization object
  non_cat_obj <- list(function_name = "cat_lmm_initialization")
  expect_false(is.cat_lmm_init(non_cat_obj))
})

test_that("is.cat_cox_init works correctly", {
  # Example of creating a cat_initialization object with correct function_name
  init_obj <- list(function_name = "cat_cox_initialization")
  class(init_obj) <- "cat_initialization"
  expect_true(is.cat_cox_init(init_obj))

  # Example with a different function_name
  init_obj_diff <- list(function_name = "cat_other_initialization")
  class(init_obj_diff) <- "cat_initialization"
  expect_false(is.cat_cox_init(init_obj_diff))

  # Example with a non-cat_initialization object
  non_cat_obj <- list(function_name = "cat_cox_initialization")
  expect_false(is.cat_cox_init(non_cat_obj))
})

test_that("is.subset_membership works correctly", {
  # Example with all elements of whole_set in the subset
  whole_set <- c(1, 2, 3)
  subset <- c(1, 2, 3)
  expect_true(is.subset_membership(whole_set, subset))

  # Example with all elements of whole_set in the subsets
  whole_set <- c(1, 2, 3)
  subset_1 <- c(1, 2)
  subset_2 <- c(3)
  expect_true(is.subset_membership(whole_set, subset_1, subset_2))

  # Example with missing elements in the subset
  whole_set <- c(1, 2, 3)
  subset <- c(1, 2)
  expect_false(is.subset_membership(whole_set, subset))

  # Example with duplicate elements in the subset
  whole_set <- c(1, 2, 3)
  subset <- c(1, 2, 2, 3)
  expect_false(is.subset_membership(whole_set, subset))

  # Example with empty subset and whole_set
  whole_set <- c()
  subset <- c()
  expect_true(is.subset_membership(whole_set, subset))

  # Example with non-empty whole_set and empty subset
  whole_set <- c(1, 2, 3)
  subset <- c()
  expect_false(is.subset_membership(whole_set, subset))
})

source("~/catalytic/tests/testthat/helper.R")

set.seed(1)

context("functions")

# ----- Functions for data_initialization -----

test_that("get_standardized_data function works correctly", {
  data <- data.frame(
    continue_col = c(1, 2, NA, 4),
    bool_col = c(1, 0, NA, 1),
    factor_col = as.factor(c("a", "b", "c", NA))
  )

  # Test case 1: Remove rows containing NA by defaut
  expect_df_equal_values(
    get_standardized_data(data),
    data.frame(
      continue_col = c(1, 2),
      bool_col = factor(c(1, 0)),
      factor_col = factor(c("a", "b"))
    )
  )

  # Test case 2: Replace NA with a specific value and mode for factor columns
  expect_df_equal_values(
    get_standardized_data(data, na_replace = 0),
    data.frame(
      continue_col = c(1, 2, 0, 4),
      bool_col = as.factor(c(1, 0, 1, 1)),
      factor_col = as.factor(c("a", "b", "c", "a"))
    )
  )

  # Test case 3: Replace NA with mean for numeric columns and mode for factor columns
  expect_df_equal_values(
    get_standardized_data(data, na_replace = stats::median),
    data.frame(
      continue_col = c(1, 2, 2, 4),
      bool_col = as.factor(c(1, 0, 1, 1)),
      factor_col = as.factor(c("a", "b", "c", "a"))
    )
  )
})

test_that("get_resampled_data performs resampling correctly", {
  # Test case 1: Balanced data
  data_balanced <- data.frame(
    col1 = c(1, 2, 3, 4, 5), # Continuous data with balanced levels
    col2 = c("A", "B", "C", "D", "E") # Categorical data with balanced levels
  )

  resampled <- get_resampled_data(data_balanced, resample_size = 10)

  expect_equal(dim(resampled$resampled_data), c(10, 2))
  expect_s3_class(resampled$resampled_data, "data.frame")
  expect_df_equal_values(
    resampled$resampled_data_record,
    data.frame(
      "Type" = c("Continuous", "Category"),
      "Process" = c("Coordinate", "Coordinate")
    )
  )

  # Test case 2: Less diversity data
  data_smooth <- data.frame(
    col1 = c(1, 2, 3, 4, 5),
    col2 = c(1, 2, 2, 3, 3) # Continuous data with three unique values and data degree two
  )

  resampled <- get_resampled_data(data_smooth,
    resample_size = 10,
    data_degree = c(1, 2)
  )
  expect_df_equal_values(
    resampled$resampled_data_record,
    data.frame(
      "Type" = c("Continuous", "Continuous"),
      "Process" = c("Coordinate", "Coordinate->Smoothing")
    )
  )

  # Test case 3: Skewed data
  data_iqr <- data.frame(
    col1 = c(1, 2, 3, 4, 5),
    col2 = c(1, 2, 2, 2, 10) # Continuous data with high skewness
  )

  resampled <- get_resampled_data(data_iqr,
    resample_size = 10
  )
  expect_df_equal_values(
    resampled$resampled_data_record,
    data.frame(
      "Type" = c("Continuous", "Continuous"),
      "Process" = c("Coordinate", "Coordinate->IQR")
    )
  )

  # Test case 4: unbalanced level data
  data_unbalance <- data.frame(
    col1 = c(1, 2, 3, 4, 5),
    col2 = c("A", "B", "B", "B", "B") # Categorical data with unbalanced levels
  )

  resampled <- get_resampled_data(data_unbalance,
    resample_size = 10
  )
  expect_df_equal_values(
    resampled$resampled_data_record,
    data.frame(
      "Type" = c("Continuous", "Category"),
      "Process" = c("Coordinate", "Coordinate->Flattening")
    )
  )
})

test_that("get_resampling_methods performs resampling correctly", {
  original_data <- c(1, rep(100, 9))


  # Test case 1: Symmetrizing method
  resampled_data <- get_resampling_methods(original_data, original_data, method = "symmetrizing")
  expect_true(length(unique(resampled_data)) > length(unique(original_data)))

  # Test case 2: Smoothing method
  resampled_data <- get_resampling_methods(original_data, original_data, method = "smoothing")
  expect_true(length(unique(resampled_data)) > length(unique(original_data)))

  # Test case 3: Flattening method
  resampled_data <- get_resampling_methods(original_data, original_data, method = "flattening")
  expect_true(sum(resampled_data == 1) > sum(original_data == 1))

  # Test case 4: IQR method
  resampled_data <- get_resampling_methods(original_data, original_data, method = "IQR")
  expect_true(abs(get_skewness(resampled_data)) < abs(get_skewness(original_data)))
})

test_that("get_adjusted_data adjusts data frame correctly", {
  data_ <- mock_data_initialization_glm_gaussian()

  # Call the function
  adjusted_data <- get_adjusted_data(data_, "y ~ x1")

  # Check if adjustments are made correctly
  expect_true("adj_x" %in% names(adjusted_data))
  expect_true("adj_obs_x" %in% names(adjusted_data))
  expect_true("adj_syn_x" %in% names(adjusted_data))

  # Check if all other part still same
  other_part_of_adjusted_data <- adjusted_data[setdiff(
    names(adjusted_data),
    c(
      "adj_x",
      "adj_obs_x",
      "adj_syn_x"
    )
  )]
  class(other_part_of_adjusted_data) <- class(data_)
  expect_equal(data_, other_part_of_adjusted_data)

  # Check adjusted data
  expect_df_equal_values(
    data.frame(x1 = data_$x[, "x1"]),
    adjusted_data$adj_x
  )
  expect_df_equal_values(
    data.frame(x1 = data_$obs_x[, "x1"]),
    adjusted_data$adj_obs_x
  )
  expect_df_equal_values(
    data.frame(x1 = data_$syn_x[, "x1"]),
    adjusted_data$adj_syn_x
  )
})


# ----- Functions for cat -----

test_that("get_stan_model retrieves correct Stan models", {
  memoise::forget(get_stan_model)

  # Test case 1: Binomial logistic regression with joint theta
  start_time <- Sys.time()
  model_glm_gaussian <- get_stan_model(model_type = "glm", glm_family = gaussian)
  time1 <- Sys.time() - start_time
  expect_s4_class(model_glm_gaussian, "stanmodel")

  # Test case 2: Memorised same input function
  start_time <- Sys.time()
  model_glm_gaussian <- get_stan_model(model_type = "glm", glm_family = gaussian)
  time2 <- Sys.time() - start_time
  expect_true(time1 > time2)

  # Test case 3: Check file names
  model_glm_gaussian_joint <- get_stan_model(
    model_type = "glm",
    glm_family = gaussian,
    joint_tau = TRUE
  )

  model_glm_binomial <- get_stan_model(
    model_type = "glm",
    glm_family = binomial
  )

  model_glm_binomial_joint <- get_stan_model(
    model_type = "glm",
    glm_family = binomial,
    joint_tau = TRUE
  )

  model_glm_binomial_joint_theta <- get_stan_model(
    model_type = "glm",
    glm_family = binomial,
    glm_binomial_joint_theta = TRUE
  )

  model_glm_binomial_joint_theta_alpha <- get_stan_model(
    model_type = "glm",
    glm_family = binomial,
    glm_binomial_joint_theta = TRUE,
    glm_binomial_joint_alpha = TRUE
  )

  model_glm_binomial_gibbs <- get_stan_model(
    model_type = "glm",
    glm_family = binomial,
    glm_binomial_gibbs = TRUE
  )


  model_cox <- get_stan_model(model_type = "cox")

  model_cox_joint <- get_stan_model(model_type = "cox", joint_tau = TRUE)


  expect_equal(attr(model_glm_gaussian@model_code, "model_name2"), "glm_gaussian")
  expect_equal(attr(model_glm_gaussian_joint@model_code, "model_name2"), "glm_gaussian_joint")
  expect_equal(attr(model_glm_binomial@model_code, "model_name2"), "glm_binomial")
  expect_equal(attr(model_glm_binomial_joint@model_code, "model_name2"), "glm_binomial_joint")
  expect_equal(attr(model_glm_binomial_joint_theta@model_code, "model_name2"), "glm_binomial_joint_theta")
  expect_equal(attr(model_glm_binomial_joint_theta_alpha@model_code, "model_name2"), "glm_binomial_joint_theta_alpha")
  expect_equal(attr(model_glm_binomial_gibbs@model_code, "model_name2"), "glm_binomial")
  expect_equal(attr(model_cox@model_code, "model_name2"), "cox")
  expect_equal(attr(model_cox_joint@model_code, "model_name2"), "cox_joint")
})

# ----- Functions for cat_lmm -----

test_that("update_sigma_square_cat preform calculation correctly", {
  sigma_square <- 0.5
  sigma_u_square <- 1.0
  d_seq <- c(1, 2, 3)
  d_seq_star <- c(1.5, 2.5, 3.5)
  q <- c(0.1, 0.2, 0.3)
  q_star <- c(0.15, 0.25, 0.35)
  tau <- 0.01

  # Test case 1: Basic calculation
  expected_result <- sum(log(sigma_u_square * d_seq + sigma_square) + q^2 / (sigma_u_square * d_seq + sigma_square)) +
    tau / length(d_seq_star) *
      sum(log(sigma_u_square * d_seq_star + sigma_square) + q_star^2 / (sigma_u_square * d_seq_star + sigma_square))

  result <- update_sigma_square_cat(sigma_square, sigma_u_square, d_seq, d_seq_star, q, q_star, tau)
  expect_equal(result, expected_result)
})

# ----- Functions for cat_cox -----

test_that("get_cox_hessian calculates Hessian matrix correctly", {
  # Example data
  X <- matrix(rnorm(50), nrow = 10)
  time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
  event <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 1)
  coefs <- rnorm(ncol(X))

  # Calculate the Hessian matrix using the function
  hessian_matrix <- get_cox_hessian(X, time, event, coefs)

  # Check that the result is a matrix of the correct dimensions
  expect_true(is.matrix(hessian_matrix))
  expect_equal(dim(hessian_matrix), c(ncol(X), ncol(X)))

  # Check that the function returns a numeric matrix
  expect_type(hessian_matrix, "double")

  # Check for specific known properties of the Hessian
  expect_true(all(diag(hessian_matrix) <= 0)) # Diagonal elements should be non-positive
})


test_that("get_syn_cox_hessian calculates synthetic Hessian matrix correctly", {
  # Example data
  X <- matrix(rnorm(50), nrow = 10)
  time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
  coefs <- rnorm(ncol(X))
  hazard_constant <- 0.1

  # Calculate the synthetic Hessian matrix using the function
  hessian_matrix <- get_syn_cox_hessian(X, time, coefs, hazard_constant)

  # Check that the result is a matrix of the correct dimensions
  expect_true(is.matrix(hessian_matrix))
  expect_equal(dim(hessian_matrix), c(ncol(X), ncol(X)))

  # Check that the function returns a numeric matrix
  expect_type(hessian_matrix, "double")

  # Check for specific known properties of the Hessian
  expect_true(all(diag(hessian_matrix) <= 0)) # Diagonal elements should be non-positive
})

test_that("get_cox_gradient calculates gradient correctly", {
  # Example data
  X <- matrix(rnorm(50), nrow = 10)
  time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
  event <- c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1)
  coefs <- rnorm(ncol(X))

  # Calculate the gradient using the function
  gradient <- get_cox_gradient(X, time, event, coefs)

  # Check that the result is a vector of the correct length
  expect_true(is.vector(gradient))
  expect_equal(length(gradient), ncol(X))

  # Check that the function returns a numeric vector
  expect_type(gradient, "double")
})

test_that("get_syn_cox_gradient calculates gradient correctly", {
  # Example data
  X <- matrix(rnorm(50), nrow = 10)
  time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
  coefs <- rnorm(ncol(X))
  hazard_constant <- 0.01

  # Calculate the gradient using the function
  gradient <- get_syn_cox_gradient(X, time, coefs, hazard_constant)

  # Check that the result is a vector of the correct length
  expect_true(is.vector(gradient))
  expect_equal(length(gradient), ncol(X))

  # Check that the function returns a numeric vector
  expect_type(gradient, "double")
})

test_that("get_cox_partial_likelihood calculates partial likelihood correctly", {
  # Example data
  X <- matrix(rnorm(50), nrow = 10)
  time <- c(5, 8, 12, 18, 22, 30, 35, 40, 42, 50)
  event <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 1)
  coefs <- rnorm(ncol(X))
  entry_points <- NULL

  # Calculate the partial likelihood using the function
  pl <- get_cox_partial_likelihood(X, time, event, coefs, entry_points)

  # Check that the result is a numeric value
  expect_type(pl, "double")
})

test_that("get_risk_and_failure_sets computes risk and failure sets correctly", {
  time_vector <- c(5, 8, 12, 18, 22)
  event_vector <- c(1, 0, 1, 1, 0)
  indicator_vector <- seq(1, 25, by = 1)

  result <- get_risk_and_failure_sets(time_vector, event_vector, indicator_vector)

  # Manually set expected values
  expected_risk_set <- matrix(0, length(time_vector), length(indicator_vector))
  expected_failure_set <- matrix(0, length(time_vector), length(indicator_vector))

  expected_risk_set[1, 1:5] <- 1
  expected_risk_set[2, 1:8] <- 1
  expected_risk_set[3, 1:12] <- 1
  expected_risk_set[4, 1:18] <- 1
  expected_risk_set[5, 1:22] <- 1

  expected_failure_set[1, 5] <- 1
  expected_failure_set[3, 12] <- 1
  expected_failure_set[4, 18] <- 1

  expect_equal(result$risk_set, expected_risk_set)
  expect_equal(result$failure_set, expected_failure_set)
})

test_that("get_risk_set_idx identifies risk set indices correctly", {
  time_of_interest <- 10
  entry_vector <- c(0, 0, 0, 5, 7)
  time_vector <- c(5, 8, 12, 18, 22)
  event_vector <- c(1, 0, 1, 1, 0)

  # Test case 1: Basic usage
  expected_risk_set_idx <- c(3, 4, 5)
  risk_set_idx <- get_risk_set_idx(time_of_interest, entry_vector, time_vector, event_vector)
  expect_equal(risk_set_idx, expected_risk_set_idx)

  # Test case 2: No subjects at risk
  time_of_interest <- 25
  expected_risk_set_idx <- integer(0)
  risk_set_idx <- get_risk_set_idx(time_of_interest, entry_vector, time_vector, event_vector)
  expect_equal(risk_set_idx, expected_risk_set_idx)

  # Test case 3: All subjects at risk
  entry_vector <- c(0, 0, 0, 0, 0)
  time_vector <- c(1, 2, 3, 4, 5)
  event_vector <- c(1, 1, 1, 0, 0)
  time_of_interest <- 1
  expected_risk_set_idx <- c(1, 2, 3, 4, 5)
  risk_set_idx <- get_risk_set_idx(time_of_interest, entry_vector, time_vector, event_vector)
  expect_equal(risk_set_idx, expected_risk_set_idx)
})

test_that("get_qr_solve function returns solution correctly", {
  A <- matrix(c(1, 1, 1, 2, 3, 5, 4, 6, 8), nrow = 3)
  b <- c(10, 20, 30)

  # Test case 1: Basic usage
  expect_equal(get_qr_solve(A, b), qr.coef(qr(A), b))

  # Test case 2: Handle NA values
  A_with_na <- matrix(c(1, 1, NA, 2, 3, 5, 4, 6, 8), nrow = 3)
  b_with_na <- c(10, 20, NA)

  # Check if there are no NA values in the solution
  expect_false(any(is.na(get_qr_solve(A_with_na, b_with_na))))
})

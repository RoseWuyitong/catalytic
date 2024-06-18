source("~/catalytic/tests/testthat/helper.R")

set.seed(1)

context("utils")

# ----- GLM Family Analysis Functions -----

test_that("validate_glm_family function validates family objects correctly", {
  # Test case 1: Validating a family object
  expect_equal(validate_glm_family(gaussian()), gaussian(link = "identity"))

  # Test case 2: Validating a family object by its name
  expect_equal(validate_glm_family("binomial"), binomial(link = "logit"))

  # Test case 3: Validating a family object using a function
  expect_equal(validate_glm_family(function() gaussian()), gaussian(link = "identity"))
})

test_that("is.binomial detects binomial model types correctly", {
  # Test case 1: Model type is "binomial"
  expect_true(is.binomial("binomial"))

  # Test case 2: Model type is not "binomial"
  expect_false(is.binomial("gaussian"))
})

test_that("is.gaussian detects gaussian model types correctly", {
  # Test case 1: Model type is "gaussian"
  expect_true(is.gaussian("gaussian"))

  # Test case 2: Model type is not "gaussian"
  expect_false(is.gaussian("binomial"))
})

test_that("get_glm_family retrieves the family name correctly", {
  # Test case 1: Using a family object
  expect_equal(get_glm_family(binomial()), "binomial")

  # Test case 2: Using a family object with link details
  expect_equal(get_glm_family(binomial(), with_link = TRUE), "binomial [logit]")

  # Test case 3: Using the name of a family
  expect_equal(get_glm_family("gaussian"), "gaussian")

  # Test case 4: Using a function that returns a family object
  expect_equal(get_glm_family(function() gaussian()), "gaussian")

  # Test case 5: Invalid family input
  expect_error(get_glm_family("invalid_family"), "'family' must be a family.")
})

# ----- GLM Mean Calculation Functions -----

test_that("mean_binomial calculates mean of binomial correctly", {
  X <- mock_data_glm_binomial()[, c("x1", "x2")]
  coefficients <- c(0.5, -0.3)
  coefficients_with_intercept <- c(0.2, 0.5, -0.3)
  X_error <- matrix(c(1, 2, 3, 4), ncol = 4)

  # Test case 1: Usage without intercept
  expect_equal(
    mean_binomial(X, coefficients),
    1 / (1 + exp(-get_linear_predictor(X, coefficients)))
  )

  # Test case 2: Usage with intercept
  expect_equal(
    mean_binomial(X, coefficients_with_intercept),
    1 / (1 + exp(-get_linear_predictor(X, coefficients_with_intercept[-1], coefficients_with_intercept[1])))
  )

  # Test case 3: Error handling with mismatched dimensions
  expect_error(
    mean_binomial(X_error, coefficients),
    "Dimensions of coefficients and predictor variables do not match."
  )
})

test_that("mean_gaussian calculates mean of gaussian correctly", {
  X <- mock_data_glm_gaussian()[, c("x1", "x2")]
  coefficients <- c(0.5, -0.3)
  coefficients_with_intercept <- c(0.2, 0.5, -0.3)
  X_error <- matrix(c(1, 2, 3, 4), ncol = 4)

  # Test case 1: Usage without intercept
  expect_equal(mean_gaussian(X, coefficients), get_linear_predictor(X, coefficients))

  # Test case 2: Usage with intercept
  expect_equal(
    mean_gaussian(X, coefficients_with_intercept),
    get_linear_predictor(X, coefficients_with_intercept[-1], coefficients_with_intercept[1])
  )

  # Test case 3: Error handling with mismatched dimensions
  expect_error(
    mean_gaussian(X_error, coefficients),
    "Dimensions of coefficients and predictor variables do not match."
  )
})

test_that("get_glm_mean calculates GLM means correctly", {
  # Generate sample data
  data <- mock_data_glm_binomial()
  X <- data[, c("x1", "x2")]
  family <- binomial

  model <- stats::glm(
    formula = y ~ x1 + x2,
    data = data,
    family = family
  )

  coefs <- coef(model)

  # Test case 1: Usage with coefficients
  glm_mean_1 <- get_glm_mean(
    X = X,
    family = family,
    coefficients = coefs
  )
  expect_equal(glm_mean_1, 1 / (1 + exp(-get_linear_predictor(
    X,
    as.vector(coefs[-1]),
    coefs[1]
  ))))

  # Test case 2: Error handling with unsupported family
  unsupported_family <- stats::poisson
  expect_error(
    get_glm_mean(X = X, family = unsupported_family),
    "Invalid GLM family."
  )
})

# ----- GLM Log Density/Log Density Gradient Calculation Functions -----

test_that("log_density_binomial calculates log density of binomail correctly", {
  data <- mock_data_glm_binomial()
  Y <- data[, "y"]
  coefs <- as.vector(coef(stats::glm(
    formula = y ~ x1 + x2,
    data = data,
    family = binomial
  )))

  lp <- get_linear_predictor(data[, c("x1", "x2")], coefs)

  weights <- c(1, 2, 1, 0.5, 2)
  weight <- 0.5

  # Test case 1: Usage with lp
  expect_equal(
    log_density_binomial(Y = Y, lp = lp),
    sum(lp * Y - log(1 + exp(lp)))
  )

  # Test case 2: Using weights for each observation
  expect_equal(
    log_density_binomial(Y = Y, lp = lp, weights = weights),
    sum(weights * (lp * Y - log(1 + exp(lp))))
  )

  # Test case 3: Using single weight for all observation
  expect_equal(
    log_density_binomial(Y, lp, weight),
    sum(weight * (lp * Y - log(1 + exp(lp))))
  )
})

test_that("log_density_gaussian calculates log density of gaussian correctly", {
  residuals <- c(1.2, -0.5, 0.7, 1.1)
  weights <- c(1, 2, 1, 0.5)
  weight <- 0.5

  # Test case 1: Usage with equal weights
  expect_equal(
    log_density_gaussian(residuals),
    0.5 * sum(1 + log(2 * pi * var(residuals)))
  )

  # Test case 2: Using weights for each residual
  expect_equal(
    log_density_gaussian(residuals, weights),
    0.5 * sum(weights * (1 + log(2 * pi * var(residuals))))
  )


  # Test case 3: Using single weight for all residuals
  expect_equal(
    log_density_gaussian(residuals, weight),
    0.5 * sum(weight * (1 + log(2 * pi * var(residuals))))
  )
})

test_that("log_density_grad_binomial calculates binomial gradient correctly", {
  X <- matrix(c(0.5, 1.2, -0.3, 0.7, 0.8, -1.1), ncol = 2)
  Y <- c(1, 0, 1)
  mean <- c(0.7, 0.4, 0.6)
  weights <- c(1, 2, 0.5)
  weight <- 0.5

  # Test case 1: Usage with equal weights
  expect_equal(
    log_density_grad_binomial(X, Y, mean),
    c(t(1 * (Y - mean)) %*% cbind(1, X))
  )

  # Test case 2: Using weights for each observation
  expect_equal(
    log_density_grad_binomial(X, Y, mean, weights),
    c(t(weights * (Y - mean)) %*% cbind(1, X))
  )

  # Test case 3: Using single weight for all observations
  expect_equal(
    log_density_grad_binomial(X, Y, mean, weight),
    c(t(weight * (Y - mean)) %*% cbind(1, X))
  )
})

# ----- Formula Analysis Functions -----

test_that("get_formula_string converts formulas to single-line strings correctly", {
  # Test case 1: Simple formula
  formula <- y ~ x1 + x2 + x3
  formula_str <- get_formula_string(formula)
  expect_equal(formula_str, "y ~ x1 + x2 + x3")

  # Test case 2: Formula with single predictor
  formula <- y ~ x1
  formula_str <- get_formula_string(formula)
  expect_equal(formula_str, "y ~ x1")

  # Test case 3: Formula with interaction term
  formula <- y ~ x1 * x2
  formula_str <- get_formula_string(formula)
  expect_equal(formula_str, "y ~ x1 * x2")

  # Test case 4: Formula with lmm formula
  formula <- y ~ x1 + (1 | group)
  formula_str <- get_formula_string(formula)
  expect_equal(formula_str, "y ~ x1 + (1 | group)")
})

test_that("get_formula_rhs extracts the right-hand side of a formula correctly", {
  # Test case 1: Extract RHS without tilde
  formula <- y ~ x1 + x2
  rhs <- get_formula_rhs(formula)
  expect_equal(rhs, "x1 + x2")

  # Test case 2: Extract RHS with tilde
  rhs_with_tilde <- get_formula_rhs(formula, with_tilde = TRUE)
  expect_equal(rhs_with_tilde, "~ x1 + x2")

  # Test case 3: Extract RHS from a formula with interaction terms
  formula <- y ~ x1 * x2
  rhs <- get_formula_rhs(formula)
  expect_equal(rhs, "x1 * x2")

  # Test case 4: Extract RHS from a formula with lmm formula
  formula <- y ~ x1 + (1 | group)
  rhs <- get_formula_rhs(formula)
  expect_equal(rhs, "x1 + (1 | group)")
})

test_that("get_formula_lhs extracts the left-hand side of a formula correctly", {
  # Test case 1: Extract LHS from a simple formula
  formula <- y ~ x1 + x2s
  lhs <- get_formula_lhs(formula)
  expect_equal(lhs, "y")

  # Test case 2: Extract LHS from a formula with a function
  formula <- log(y) ~ x1 + x2
  lhs <- get_formula_lhs(formula)
  expect_equal(lhs, "log(y)")

  # Test case 3: Extract LHS from a formula with nested functions
  formula <- sqrt(log(y)) ~ x1 + x2
  lhs <- get_formula_lhs(formula)
  expect_equal(lhs, "sqrt(log(y))")

  # Test case 4: Extract LHS from a formula with an interaction term
  formula <- y1:y2 ~ x1 + x2
  lhs <- get_formula_lhs(formula)
  expect_equal(lhs, "y1:y2")
})

test_that("is.formula detects formulas correctly", {
  # Test case 1: Check if a formula object is detected correctly
  expect_true(is.formula(as.formula("y ~ x1 + x2")))

  # Test case 2: Check if a survival formula object is detected correctly
  expect_true(is.formula(as.formula("Surv(time, status) ~ 1")))

  # Test case 3: Check if a lmm formula object is detected correctly
  expect_true(is.formula(as.formula("y ~ 1 + (x|z)")))

  # Test case 4: Check if a lmm formula object is detected correctly
  expect_true(is.formula(as.formula("~ 1 + x1")))

  # Test case 5: Check if a non-formula object is detected correctly
  expect_false(is.formula(c(1, 2, 3)))
})

test_that("is.glm_formula detects GLM formulas correctly", {
  # Test case 1: GLM formula with predictors
  expect_true(is.glm_formula(as.formula("y ~ x1 + x2 + x3")))

  # Test case 2: GLM formula with only an intercept
  expect_true(is.glm_formula(as.formula("y ~ 1")))

  # Test case 3: Non-GLM formula
  expect_false(is.glm_formula(as.formula("Surv(y, z)~ x1 + x2")))

  # Test case 4: Non-GLM formula (no response)
  expect_false(is.glm_formula(as.formula("~ x1 + x2")))
})

test_that("is.lmm_formula detects LMM formulas correctly", {
  # Test case 1: LMM formula with fixed and random effects
  expect_true(is.lmm_formula(as.formula("y ~ x1 + (1 | subject)")))

  # Test case 2: Non-LMM formula
  expect_false(is.lmm_formula(as.formula("y ~ x1 + x2")))
})

test_that("is.surv_formula detects survival formulas correctly", {
  # Test case 1: Formula for survival analysis
  expect_true(is.surv_formula(as.formula("Surv(time, status) ~ age + treatment")))

  # Test case 2: Non-survival formula
  expect_false(is.surv_formula(as.formula("y ~ x1 + x2")))
})

# ----- Data Initialization Type/Class Identification Functions -----

data_glm_gaussian <- mock_data_initialization_glm_gaussian()
data_glm_binomial <- mock_data_initialization_glm_binomial()
data_lmm <- mock_data_initialization_lmm()
data_surv <- mock_data_initialization_surv()

test_that("is.data_glm detects data generated from data_initilization_glm correctly", {
  # Test case 1: Data object is generated from data_initilization_glm
  expect_true(is.data_glm(data_glm_gaussian))
  expect_true(is.data_glm(data_glm_binomial))

  # Test case 2: Data object is not generated from data_initilization_glm
  expect_false(is.data_glm(data_lmm))
})

test_that("is.data_lmm detects data generated from data_initilization_lmm correctly", {
  # Test case 1: Data object is generated from data_initilization_glm
  expect_true(is.data_lmm(data_lmm))

  # Test case 2: Data object is not generated from data_initilization_glm
  expect_false(is.data_lmm(data_surv))
})

test_that("is.data_surv detects data generated from data_initilization_surv correctly", {
  # Test case 1: Data object is generated from data_initilization_glm
  expect_true(is.data_surv(data_surv))

  # Test case 2: Data object is not generated from data_initilization_glm
  expect_false(is.data_surv(data_glm_gaussian))
})

# ----- Discrepancy Functions -----

test_that("discrepancy_logarithmic_error calculates logarithmic error correctly", {
  Y <- c(1.3, 2, 3, 4)
  est_Y <- c(1.1, 2.0, 2.8, 4.2)

  # Test case 1: Basic usage
  expect_equal(
    discrepancy_logarithmic_error(Y, est_Y),
    Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y))
  )
})

test_that("discrepancy_square_error calculates square error correctly", {
  Y <- c(1, 2, 3, 4)
  est_Y <- c(0.9, 2.1, 2.8, 4.2)

  # Test case 1: Basic usage
  expect_equal(discrepancy_square_error(Y, est_Y), (Y - est_Y)^2)
})

test_that("discrepancy_classification_error calculates classification error correctly", {
  Y <- c(1, 0, 1, 1)
  est_Y <- c(1, 0, 1, 0)

  # Test case 1: Basic usage
  expect_equal(
    discrepancy_classification_error(Y, est_Y),
    as.numeric(Y == est_Y) / length(Y)
  )
})

test_that("get_discrepancy calculates discrepancy correctly for various methods", {
  binomial_Y <- c(1, 0, 1, 1)
  binomial_response_est_Y <- c(1, 0, 1, 1)
  binomial_prob_est_Y <- c(0.1, 0.3, 0.2, 0.5)

  gaussian_Y <- c(2.2, 3.3, 4.4, 5.5, 6)
  gaussian_est_Y <- c(2, 3, 4, 5, 6)

  # Test case 1: Classification error
  expect_equal(
    get_discrepancy(
      discrepancy_method = "classification_error",
      Y = binomial_Y,
      est_Y = binomial_response_est_Y,
      family = binomial
    ),
    as.numeric(binomial_Y == binomial_response_est_Y) / length(binomial_Y)
  )

  # Test case 2: Classification error - error handling with wrong data input
  expect_error(
    get_discrepancy(
      discrepancy_method = "classification_error",
      Y = gaussian_Y,
      est_Y = gaussian_est_Y,
      family = binomial
    ),
    "`classification_error` can only use for 0 or 1 input."
  )

  # Test case 3: Classification error - error handling wrong family input
  expect_error(
    get_discrepancy(
      discrepancy_method = "classification_error",
      Y = Y,
      est_Y = est_Y,
      family = gaussian
    ),
    "`classification_error` can only use for `binomial` family"
  )


  # Test case 4: Logarithmic error
  expect_equal(
    get_discrepancy(
      discrepancy_method = "logarithmic_error",
      Y = gaussian_Y,
      est_Y = gaussian_est_Y
    ),
    discrepancy_logarithmic_error(
      Y = pmax(0.0001, pmin(0.9999, gaussian_Y)),
      est_Y = pmax(0.0001, pmin(0.9999, gaussian_est_Y))
    )
  )

  # Test case 5: Square error
  expect_equal(
    get_discrepancy(
      discrepancy_method = "square_error",
      Y = gaussian_Y,
      est_Y = gaussian_est_Y
    ),
    (gaussian_Y - gaussian_est_Y)^2
  )

  # Test case 6: Square error - warning handling wrong family input
  expect_warning(
    get_discrepancy(
      discrepancy_method = "square_error",
      Y = gaussian_Y,
      est_Y = gaussian_est_Y,
      family = binomial
    ),
    "Recommand use `square_error` only `gaussian` family"
  )

  # Test case 7: Logistic deviance
  expect_equal(
    get_discrepancy(
      discrepancy_method = "logistic_deviance",
      Y = binomial_Y,
      est_Y = binomial_prob_est_Y,
      family = binomial
    ),
    -log_density_binomial(Y = binomial_Y, lp = log(1 / binomial_prob_est_Y - 1))
  )

  # Test case 8: Logistic deviance - error handling with wrong data input
  expect_error(
    get_discrepancy(
      discrepancy_method = "logistic_deviance",
      Y = binomial_Y,
      est_Y = binomial_response_est_Y,
      family = binomial
    ),
    "`logistic_deviance` can only use when estimated Y is probability."
  )

  # Test case 9: Logistic deviance - error handling with wrong family input
  expect_error(
    get_discrepancy(
      discrepancy_method = "logistic_deviance",
      Y = binomial_Y,
      est_Y = binomial_prob_est_Y,
      family = gaussian
    ),
    "`logistic_deviance` can only use for `binomial` family"
  )
})

# ----- Miscellaneous Functions -----

test_that("is.continuous detects continuous variables correctly", {
  # Test case 1: Numeric vector with more than two unique values
  expect_true(is.continuous(c(1.5, 2.3, 3.7, 4.2)))

  # Test case 2: Numeric vector with only two unique values
  expect_false(is.continuous(c(1.5, 2.5, 1.5, 2.5)))

  # Test case 3: Non-numeric vector
  expect_false(is.continuous(c("a", "b", "c")))
})

test_that("get_linear_predictor calculates linear predictor correctly", {
  X <- matrix(c(1, 2, 3, 4), ncol = 2)
  coefs <- c(0.5, -0.3)
  coefs_with_intercept <- c(0.2, 0.5, -0.3)

  # Test case 1: Usage without intercept
  expect_equal(get_linear_predictor(X, coefs), X %*% coefs)

  # Test case 2: Usage with intercept
  expect_equal(
    get_linear_predictor(X, coefs_with_intercept),
    X %*% coefs_with_intercept[-1] + coefs_with_intercept[1]
  )
})

# test_that("get_correlation_matrix constructs correlation matrix correctly", {
#   correlation_coeffs <- c(0.8, 0.3, 0.5)
#   dimension <- 3
#   expected_corr_matrix <- matrix(
#     c(1.0, 0.3, 0.5,
#       0.8, 1.0, 0.5,
#       0.8, 0.3, 1.0),
#     nrow = 3, ncol = 3)
#   actual_corr_matrix <- get_correlation_matrix(correlation_coeffs, dimension)
#   expect_equal(actual_corr_matrix, expected_corr_matrix)
# })

test_that("check_subset_membership verifies subset membership correctly", {
  # Test case 1: Valid subset membership
  expect_true(check_subset_membership(1:5, 1:3, 4:5))

  # Test case 2: Subset overlap
  expect_false(check_subset_membership(1:5, 1:3, 3:5))

  # Test case 3: Missing elements in subset
  expect_false(check_subset_membership(1:5, 1:4))
})

test_that("get_iqr function computes IQR correctly", {
  x <- c(10, 20, 30, 40, 50)
  # Test case 1: Basic usage with default interval
  expect_equal(get_iqr(x), 20)

  # Test case 2: Usage with different intervals
  expect_equal(get_iqr(x, 0.25), 30)

  # Test case 3: Edge case with single value
  x_single <- c(10)
  expected_iqr_single <- 0
  expect_equal(get_iqr(x_single), 0)

  # Test case 4: Edge case with repeated values
  x_repeat <- c(10, 10, 10, 10, 50)
  expect_equal(get_iqr(x_repeat), 8)
})

test_that("get_lambda_function returns correct lambda functions", {
  X <- matrix(rnorm(100), ncol = 5)
  coefficients <- rnorm(5)

  # Test case 1: square_error
  lambda_square_error <- get_lambda_function("square_error", X, coefficients)
  expect_equal(lambda_square_error, 2 * get_glm_mean(family = "gaussian", X = X, coefficients = coefficients))

  # Test case 2: classification_error
  lambda_classification_error <- get_lambda_function("classification_error", X, coefficients)
  expect_equal(lambda_classification_error, 2 * get_glm_mean(family = "binomial", X = X, coefficients = coefficients))

  # Test case 3: logistic_deviance
  lambda_logistic_deviance <- get_lambda_function("logistic_deviance", X, coefficients)
  expect_equal(lambda_logistic_deviance, get_linear_predictor(X, coefficients))
})

test_that("get_mode function works correctly", {
  # Test case 1: Mode of a numeric vector
  numeric_vector <- c(1, 2, 2, 3, 4)
  expect_equal(get_mode(numeric_vector), 2)

  # Test case 2: Mode of a character vector
  char_vector <- c("a", "b", "a", "c")
  expect_equal(get_mode(char_vector), "a")

  # Test case 3: Mode of a factor vector
  factor_vector <- factor(c("low", "medium", "medium", "high"))
  expect_equal(as.character(get_mode(factor_vector)), "medium")

  # Test case 4: Vectors with a single unique value
  single_value_vector <- c(5, 5, 5, 5)
  expect_equal(get_mode(single_value_vector), 5)

  # Test case 5: Vectors with NA values
  na_vector <- c(1, 2, 2, 3, NA)
  expect_equal(get_mode(na_vector), 2)

  # Test case 6: Returns the first mode in case of a tie
  tie_vector <- c(1, 2, 2, 3, 3)
  expect_equal(get_mode(tie_vector), 2)
})

test_that("get_skewness calculates skewness correctly", {
  vec1 <- c(1, 2, 3, 2, 1)
  vec2 <- c(1, 1, 1, 1, 1)

  expect_equal(get_skewness(vec1), sum((vec1 - mean(vec1))^3) / (length(vec1) * sd(vec1)^3))
  expect_equal(get_skewness(vec2), NaN)
})

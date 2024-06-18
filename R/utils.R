# ----- GLM Family Analysis Functions -----

#' Validates the family object for GLM
#'
#' This function validates the family object for use in a generalized linear model (GLM).
#' It ensures that the provided family object is of class "family".
#'
#' @param f A character string representing the name of the family object or
#'          a function that returns a family object.
#'
#' @return A valid family object for GLM.
#' @export
#'
#' @examples
#' # Example 1: Validating a family object
#' my_family <- gaussian()
#' validated_family <- validate_glm_family(my_family)
#' print(validated_family)
#'
#' # Example 2: Validating a family object by its name
#' family_name <- "binomial"
#' validated_family <- validate_glm_family(family_name)
#' print(validated_family)
#'
#' # Example 3: Validating a family object using a function
#' get_family <- function() gaussian()
#' validated_family <- validate_glm_family(get_family)
#' print(validated_family)
validate_glm_family <- function(f) {
  if (is.character(f)) {
    tryCatch(
      {
        f <- get(f, mode = "function", envir = parent.frame(2))
      },
      error = function(.) stop("'family' must be a family.", call. = FALSE)
    )
  }
  if (is.function(f)) {
    f <- f()
  }
  if (!is(f, "family")) {
    stop("'family' must be a family.", call. = FALSE)
  }

  return(f)
}

#' Checks if the model type is binomial
#'
#' This function checks if the given model type is "binomial".
#'
#' @param x The model type to be checked.
#'
#' @return TRUE if the model type is "binomial", otherwise FALSE.
#' @export
#'
#' @examples
#' is_binomial <- is.binomial("binomial")
#' print(is_binomial) # Should return TRUE
is.binomial <- function(x) {
  return(x == "binomial")
}

#' Checks if the model type is Gaussian
#'
#' This function checks if the given model type is "gaussian".
#'
#' @param x The model type to be checked.
#'
#' @return TRUE if the model type is "gaussian", otherwise FALSE.
#' @export
#'
#' @examples
#' is_gaussian <- is.gaussian("gaussian")
#' print(is_gaussian) # Should return TRUE
is.gaussian <- function(x) {
  return(x == "gaussian")
}

#' Retrieve the GLM family name with optional link function details
#'
#' This function validates the provided GLM family object and retrieves its family name.
#' Optionally, it appends details about the link function to the family name if `with_link = TRUE`.
#'
#' @param family A character string representing the name of the family object,
#'               a function that returns a family object, or a family object itself.
#' @param with_link Logical indicating whether to include details about the link function (`TRUE`) or not (`FALSE`).
#'                  Default is FALSE.
#'
#' @return A character string representing the GLM family name, with optional link function details if requested.
#' @export
#'
#' @examples
#' # Example 1: Retrieving the family name without link details
#' fam <- get_glm_family("binomial")
#' print(fam) # Should print "binomial"
#'
#' # Example 2: Retrieving the family name with link details
#' fam <- get_glm_family("binomial", with_link = TRUE)
#' print(fam) # Should print "binomial [logit]"
#'
#' # Example 3: Using a function that returns a family object
#' fam <- get_glm_family(function() binomial())
#' print(fam) # Should print "binomial"
#'
#' # Example 4: Using a function that returns a family object
#' fam <- get_glm_family(binomial(link = logit))
#' print(fam) # Should print "binomial"
get_glm_family <- function(family, with_link = FALSE) {
  family <- validate_glm_family(family)
  if (with_link) {
    return(paste0(family$family, " [", family$link, "]"))
  }
  return(family$family)
}

# ----- GLM Mean Calculation Functions -----

#' Calculates the Mean of a binomial distribution
#'
#' This function calculates the mean of a binomial distribution given the predictor variables (X) and coefficients.
#' The coefficients can include an intercept term, which is handled accordingly.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param coefficients A numeric vector of coefficients corresponding to each predictor variable.
#'                     If an intercept term is included, it should be the first element of the vector.
#'
#' @return A numeric vector representing the mean of the binomial distribution for each observation.
#' @export
#'
#' @examples
#' # Example 1: Basic usage without intercept
#' X <- matrix(c(1, 2, 3, 4), ncol = 2)
#' coefficients <- c(0.5, -0.3)
#' mean_binomial <- mean_binomial(X, coefficients)
#' print(mean_binomial)
#'
#' # Example 2: Usage with intercept
#' coefficients_with_intercept <- c(0.2, 0.5, -0.3)
#' mean_binomial_with_intercept <- mean_binomial(X, coefficients_with_intercept)
#' print(mean_binomial_with_intercept)
#'
#' # Example 3: Error handling with mismatched dimensions
#' X_error <- matrix(c(1, 2, 3, 4), ncol = 1)
#' coefficients_error <- c(0.2, 0.5, -0.3)
#' mean_binomial_error <- mean_binomial(X_error, coefficients_error)
#' print(mean_binomial_error)
mean_binomial <- function(X,
                          coefficients) {
  return(1 / (1 + exp(-get_linear_predictor(as.data.frame(X), as.vector(coefficients)))))
}

#' Calculates the Mean of a Gaussian distribution
#'
#' This function calculates the mean of a Gaussian (normal) distribution given the predictor variables (X) and coefficients.
#' The coefficients can include an intercept term, which is handled accordingly.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param coefficients A numeric vector of coefficients corresponding to each predictor variable.
#'                     If an intercept term is included, it should be the first element of the vector.
#'
#' @return A numeric vector representing the mean of the Gaussian distribution for each observation.
#' @export
#'
#' @examples
#' # Example 1: Basic usage without intercept
#' X <- matrix(c(1, 2, 3, 4), ncol = 2)
#' coefficients <- c(0.5, -0.3)
#' mean_gaussian <- mean_gaussian(X, coefficients)
#' print(mean_gaussian)
#'
#' # Example 2: Usage with intercept
#' coefficients_with_intercept <- c(0.2, 0.5, -0.3)
#' mean_gaussian_with_intercept <- mean_gaussian(X, coefficients_with_intercept)
#' print(mean_gaussian_with_intercept)
#'
#' # Example 3: Error handling with mismatched dimensions
#' X_error <- matrix(c(1, 2, 3, 4), ncol = 1)
#' coefficients_error <- c(0.2, 0.5, -0.3)
#' mean_gaussian_error <- mean_gaussian(X_error, coefficients_error)
#' print(mean_gaussian_error)
mean_gaussian <- function(X, # TODO:: change to prediction related
                          coefficients) {
  return(get_linear_predictor(as.data.frame(X), as.vector(coefficients)))
}

#' Get Mean Function for GLM Families
#'
#' This function calculates the mean function for various families in generalized linear models (GLMs).
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param family The family of the GLM model. It can be either binomial or gaussian.
#' @param model An optional GLM model object. Default is NULL.
#' @param coefficients An optional numeric vector of coefficients. Default is NULL.
#'                     If provided, these coefficients will be used instead of extracting them from the model object.
#'
#' @return A numeric vector representing the mean function for each observation.
#' @export
#'
#' @examples
#' data <- as.data.frame(matrix(
#'   data = rnorm(15),
#'   ncol = 3,
#'   dimnames = list(c(), c("y", "x1", "x2"))
#' ))
#'
#' model <- glm(y ~ x1 + x2, data = data, family = gaussian)
#' coefficients <- coef(model)
#'
#' get_glm_mean(data[, c("x1", "x2")], gaussian, coefficients = coefficients)
get_glm_mean <- function(family,
                         X,
                         coefficients = NULL) {
  switch(get_glm_family(family),
    "binomial" = function(X,
                          coefficients) {
      return(
        mean_binomial(
          X = X,
          coefficients = coefficients
        )
      )
    },
    "gaussian" = function(X,
                          coefficients) {
      return(
        mean_gaussian(
          X = X,
          coefficients = coefficients
        )
      )
    },
    stop("Invalid GLM family. Supported families are 'binomial' and 'gaussian'.", call. = FALSE)
  )(X, coefficients)
}

# ----- GLM Log Density/Log Density Gradient Calculation Functions -----

#' Calculates the log density of the binomial distribution
#'
#' This function calculates the log density of the binomial distribution based on
#' the observed outcomes (Y), lp, and optional weights.
#'
#' @param Y A numeric vector of observed outcomes, usually 0 or 1.
#' @param lp A numeric vector representing the linear predictor.
#' @param weights A numeric vector of optional weights associated with each observation. Default is 1.
#'
#' @return The log density of the binomial distribution.
#' @export
#'
#' @examples
#' # Example 1: Calculating log density using linear predictor
#' Y <- c(0, 1, 1, 0)
#' lp <- c(0.2, 0.8, 0.6, 0.3)
#' log_density <- log_density_binomial(Y, lp = lp)
#' print(log_density)
#'
#' # Example 2: Calculating log density with weight vector
#' weights <- c(1, 2, 1, 3)
#' log_density <- log_density_binomial(Y, lp = lp, weights = weights)
#' print(log_density)
#'
#' # Example 3: Calculating log density with single weight
#' weights <- 2
#' log_density <- log_density_binomial(Y, lp = lp, weights = weights)
#' print(log_density)
log_density_binomial <- function(Y,
                                 lp,
                                 weights = 1) {
  return(sum(weights * (lp * Y - log(1 + exp(lp)))))
}


#' Calculates the log density of a Gaussian distribution for
#' a given set of residuals with weight(s)
#'
#' This function is useful in the context of statistical modeling,
#' particularly in linear models where the Gaussian distribution is assumed.
#'
#' @param residuals A numeric vector of residuals,
#'                  typically the differences between observed and predicted values.
#' @param weights An optional numeric vector of weights for each observation or
#'                just a single number for all observation. The default value is 1,
#'                meaning each observation is equally weighted.
#'
#' @return A single numeric value representing the summed log density of
#'         the Gaussian distribution for the provided residuals.
#' @export
#'
#' @examples
#' # Example 1: Basic usage with equal weights
#' residuals <- c(1.2, -0.5, 0.7, 1.1)
#' log_density <- log_density_gaussian(residuals)
#' print(log_density)
#'
#' # Example 2: Using weights for each residual
#' weights <- c(1, 2, 1, 0.5)
#' log_density_weighted_vector <- log_density_gaussian(residuals, weights)
#' print(log_density_weighted_vector)
#'
#' # Example 3: Using single weight for all residuals
#' weight <- 0.5
#' log_density_weighted_single <- log_density_gaussian(residuals, weight)
#' print(log_density_weighted_single)
log_density_gaussian <- function(residuals,
                                 weights = 1) {
  return(0.5 * sum(weights * (1 + log(2 * pi * var(residuals)))))
}


#' Calculates the gradient of the log density for a binomial distribution
#' given the predictor variables, observed outcomes, predicted mean values, and weights.
#'
#' This function is useful in the context of statistical modeling, particularly
#' in generalized linear models (GLMs) where
#' the binomial distribution is used, such as logistic regression.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param Y A numeric vector of observed binary outcomes. Each element should be 0 or 1,
#'          representing the two possible outcomes of a binomial experiment.
#' @param lp A numeric vector of the same length as Y,
#'             representing the predicted lp values.
#' @param weights An optional numeric vector of weights for each observation or
#'                just a single number for all observation. The default value is 1,
#'                meaning each observation is equally weighted.
#'
#' @return A numeric vector representing the gradient of the log density
#'         with respect to the model parameters.
#' @export
#'
#' @examples
#' # Example 1: Basic usage with equal weights
#' X <- matrix(c(0.5, 1.2, -0.3, 0.7, 0.8, -1.1), ncol = 2)
#' Y <- c(1, 0, 1)
#' lp <- c(0.7, 0.4, 0.6)
#' gradient <- log_density_grad_binomial(X, Y, lp)
#' print(gradient)
#'
#' # Example 2: Using weights for each observation
#' weights <- c(1, 2, 0.5)
#' gradient_weighted_vector <- log_density_grad_binomial(X, Y, lp, weights)
#' print(gradient_weighted_vector)
#'
#' # Example 3: Using single weight for all observations
#' weight <- 0.5
#' gradient_weighted_single <- log_density_grad_binomial(X, Y, lp, weight)
#' print(gradient_weighted_single)
log_density_grad_binomial <- function(X,
                                      Y,
                                      lp,
                                      weights = 1) {
  return(c(t(weights * (Y - lp)) %*% cbind(1, X)))
}

# ----- Formula Analysis Functions -----

#' Converts a formula to a single-line string
#'
#' This function converts a given formula to a single-line string without any line breaks or extra spaces.
#'
#' @param formula A formula object to be converted to a string.
#'
#' @return A character string representing the formula.
#' @export
#'
#' @examples
#' # Example 1: Simple formula
#' formula_str <- get_formula_string(y ~ x1 + x2 + x3)
#' print(formula_str) # Should print "y ~ x1 + x2 + x3"
#'
#' # Example 2: Formula with single predictor
#' formula_str <- get_formula_string(y ~ x1)
#' print(formula_str) # Should print "y ~ x1"
get_formula_string <- function(formula) {
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  return(char)
}

#' Extracts the right-hand side of a formula
#'
#' This function extracts the right-hand side of a given formula. Optionally, it can include the tilde (~) in the output.
#'
#' @param formula A formula object from which to extract the right-hand side.
#' @param with_tilde A logical value indicating whether to include the tilde (~) in the output. Default is FALSE.
#'
#' @return If `with_tilde` is TRUE, returns a formula object representing the right-hand side with the tilde (~).
#'         If `with_tilde` is FALSE, returns a character string representing the right-hand side without the tilde (~).
#' @export
#'
#' @examples
#' # Example 1: Extract RHS without tilde
#' formula <- y ~ x1 + x2
#' rhs <- get_formula_rhs(formula)
#' print(rhs) # Should print "x1 + x2"
#'
#' # Example 2: Extract RHS with tilde
#' rhs_with_tilde <- get_formula_rhs(formula, with_tilde = TRUE)
#' print(rhs_with_tilde) # Should print "~ x1 + x2"
get_formula_rhs <- function(formula, with_tilde = FALSE) {
  if (with_tilde) {
    return(paste("~", deparse(rlang::f_rhs(formula)), collapse = " "))
  }
  return(get_formula_string(rlang::f_rhs(formula)))
}

#' Extracts the left-hand side of a formula as a string
#'
#' This function extracts the left-hand side of a given formula and converts it to a single-line string.
#'
#' @param formula A formula object from which to extract the left-hand side.
#'
#' @return A character string representing the left-hand side of the formula.
#' @export
#'
#' @examples
#' # Example 1: Extract LHS from a formula
#' formula <- y ~ x1 + x2
#' lhs <- get_formula_lhs(formula)
#' print(lhs) # Should print "y"
#'
#' # Example 2: Extract LHS from a more complex formula
#' formula <- log(y) ~ x1 + x2
#' lhs <- get_formula_lhs(formula)
#' print(lhs) # Should print "log(y)"
get_formula_lhs <- function(formula) {
  return(get_formula_string(rlang::f_lhs(formula)))
}


#' Checks if an object is a formula
#'
#' This function checks whether the input object is a formula by searching for the tilde (~) character,
#' which is typically present in R formulas.
#'
#' @param f The object to be checked.
#' @return A logical value indicating whether the object is a formula.
#' @export
#' @examples
#' # Example 1: Check if a formula object is detected correctly
#' formula_obj <- as.formula("y ~ x1 + x2")
#' is.formula(formula_obj) # Should return TRUE
#'
#' # Example 2: Check if a half-formula object is detected correctly
#' formula_obj <- as.formula("~ x1 + x2")
#' is.formula(formula_obj) # Should return TRUE
#'
#' # Example 3: Check if a non-formula object is detected correctly
#' non_formula_obj <- c(1, 2, 3)
#' is.formula(non_formula_obj) # Should return FALSE
is.formula <- function(f) {
  return(grepl("~", deparse(f)))
}


#' Checks if a formula is for a generalized linear model (GLM)
#'
#' This function checks if a given formula is suitable for a generalized linear model (GLM)
#' by ensuring it contains a mean and predictors, which is the typical structure of a GLM formula.
#'
#' @param f A formula object to be checked.
#'
#' @return A logical value indicating whether the formula is suitable for a GLM.
#' @export
#'
#' @examples
#' # Example 1: GLM formula with predictors
#' formula_glm <- as.formula("y ~ x1 + x2 + x3")
#' is_glm <- is.glm_formula(formula_glm)
#' print(is_glm) # Should return TRUE
#'
#' # Example 2: GLM formula with only an intercept
#' formula_intercept <- as.formula("y ~ 1")
#' is_glm <- is.glm_formula(formula_intercept)
#' print(is_glm) # Should return TRUE
#'
#' # Example 3: Non-GLM formula (no response)
#' formula_non_glm <- as.formula("Surv(y, z)~ x1 + x2")
#' is_glm <- is.glm_formula(formula_non_glm)
#' print(is_glm) # Should return FALSE
is.glm_formula <- function(f) {
  if (is.lmm_formula(f) | is.surv_formula(f)) {
    return(FALSE)
  }
  formula_parts <- strsplit(deparse(f), "~")[[1]]
  return(length(formula_parts) == 2 && nchar(trimws(formula_parts[1])) > 0)
}

#' Checks if a formula is for a linear mixed-effects model (LMM)
#'
#' This function checks if a given formula is suitable for a linear mixed-effects model (LMM)
#' by ensuring it contains both fixed and random effects specifications, which is the typical structure of an LMM formula.
#'
#' @param f A formula object to be checked.
#'
#' @return A logical value indicating whether the formula is suitable for an LMM.
#' @export
#'
#' @examples
#' # Example 1: LMM formula with fixed and random effects
#' formula_lmm <- as.formula("y ~ x1 + (1 | subject)")
#' is_lmm <- is.lmm_formula(formula_lmm)
#' print(is_lmm) # Should return TRUE
#'
#' # Example 2: Non-LMM formula (no random effects)
#' formula_non_lmm <- as.formula("y ~ x1 + x2")
#' is_lmm <- is.lmm_formula(formula_non_lmm)
#' print(is_lmm) # Should return FALSE
is.lmm_formula <- function(f) {
  return(grepl("\\(.*|.*\\)", deparse(f)))
}

#' Checks if a formula is for a survival analysis model
#'
#' This function checks if a given formula is suitable for survival analysis by looking for
#' the presence of terms typically associated with survival analysis, such as the "Surv" function.
#'
#' @param f A formula object to be checked.
#'
#' @return A logical value indicating whether the formula is suitable for survival analysis.
#' @export
#'
#' @examples
#' # Example 1: Formula for survival analysis
#' formula_survival <- as.formula("Surv(time, status) ~ age + treatment")
#' is_surv_formula <- is.surv_formula(formula_survival)
#' print(is_surv_formula) # Should return TRUE
#'
#' # Example 2: Non-survival formula
#' formula_non_survival <- as.formula("y ~ x1 + x2")
#' is_surv_formula <- is.surv_formula(formula_non_survival)
#' print(is_surv_formula) # Should return FALSE
is.surv_formula <- function(f) {
  return(grepl("Surv", deparse(f)))
}

# ----- Data Initialization Type/Class Identification Functions -----

#' Checks if the data object is generated from 'data_initilization_glm'
#'
#' This function checks if the given data object is generated from' data_initilization_glm'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_glm', otherwise FALSE.
#' @export
is.data_glm <- function(x) {
  return(inherits(x, "data_initilization") & x$function_name == "data_initilization_glm")
}

#' Checks if the data object is generated from 'data_initilization_lmm'
#'
#' This function checks if the given data object is generated from' data_initilization_lmm'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_lmm', otherwise FALSE.
#' @export
is.data_lmm <- function(x) {
  return(inherits(x, "data_initilization") & x$function_name == "data_initilization_lmm")
}

#' Checks if the data object is generated from 'data_initilization_cox'
#'
#' This function checks if the given data object is generated from' data_initilization_cox'.
#'
#' @param x The data object to be checked.
#'
#' @return TRUE if the data is generated from 'data_initilization_cox', otherwise FALSE.
#' @export
is.data_cox <- function(x) {
  return(inherits(x, "data_initilization") & x$function_name == "data_initilization_cox")
}

# ----- Discrepancy Functions -----

#' Calculates the discrepancy of logarithmic error
#'
#' This function computes the logarithmic error between observed values and estimated values.
#' It is useful for evaluating the performance of predictive models, particularly in contexts where
#' the logarithmic loss (log-loss) is a suitable metric.
#'
#' @param Y A numeric vector of observed values. These values should be greater than 0 to avoid logarithm of zero issues.
#' @param est_Y A numeric vector of the same length as Y representing the estimated values.
#'
#' @return A numeric vector representing the logarithmic error for each observation.
#' @export
#'
#' @examples
#' Y <- c(1.3, 2, 3, 4)
#' est_Y <- c(1.1, 2.0, 2.8, 4.2)
#' log_error <- discrepancy_logarithmic_error(Y, est_Y)
#' print(log_error)
discrepancy_logarithmic_error <- function(Y, est_Y) {
  return(Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y)))
}

#' Calculates the discrepancy of squared error
#'
#' This function computes the squared error between observed values and estimated values.
#' It is useful for evaluating the performance of predictive models, particularly in contexts where
#' the mean squared error (MSE) or similar metrics are suitable measures of model accuracy.
#'
#' @param Y A numeric vector of observed values.
#' @param est_Y A numeric vector of the same length as Y representing the estimated values.
#'
#' @return A numeric vector representing the squared error for each observation.
#' @export
#'
#' @examples
#' Y <- c(1, 2, 3, 4)
#' est_Y <- c(0.9, 2.1, 2.8, 4.2)
#' squared_error <- discrepancy_square_error(Y, est_Y)
#' print(squared_error)
discrepancy_square_error <- function(Y, est_Y) {
  return((Y - est_Y)^2)
}

#' Calculates the discrepancy of classification error
#'
#' This function computes the classification error between observed binary outcomes and estimated binary outcomes.
#' It is useful for evaluating the performance of predictive models in classification tasks, providing the proportion
#' of correctly classified instances.
#'
#' @param Y A numeric vector of observed binary outcomes. Each element should be 0 or 1,
#'          representing the two possible outcomes of a classification task.
#' @param est_Y A numeric vector of the same length as Y representing the estimated binary outcomes.
#'
#' @return A numeric value representing the proportion of correctly classified instances.
#' @export
#'
#' @examples
#' Y <- c(1, 0, 1, 1)
#' est_Y <- c(1, 0, 1, 0)
#' classification_error <- discrepancy_classification_error(Y, est_Y)
#' print(classification_error)
discrepancy_classification_error <- function(Y, est_Y) {
  return(as.numeric(Y == est_Y) / length(Y))
}


#' Calculates the discrepancy between observed and estimated values
#'
#' This function calculates the discrepancy between observed values (Y) and estimated values (est_Y) using various methods such as logarithmic error, square error, classification error, or logistic deviance. The method is selected using the `discrepancy_method` parameter.
#'
#' @param discrepancy_method A character vector specifying the method to calculate the discrepancy. Possible values are "logarithmic_error", "square_error", "classification_error", and "logistic_deviance".
#' @param Y A numeric vector of observed values.
#' @param est_Y A numeric vector of estimated values. If NULL, estimated values are computed using the provided `family`, `X`, `est_model`, and `est_coefficients`.
#' @param family The GLM family used for estimation. Default is the gaussian family.
#' @param X A matrix or data frame of predictor variables used to compute estimated values if `est_Y` is NULL.
#' @param est_model The estimated model object used to compute estimated values if `est_Y` is NULL.
#' @param est_coefficients The estimated coefficients used to compute estimated values if `est_Y` is NULL.
#'
#' @return The discrepancy between observed and estimated values based on the selected method.
#' @export
#'
#' @examples
#' Y <- c(1, 0, 1, 1)
#' est_Y <- c(1, 0, 1, 0)
#' discrepancy <- get_discrepancy(discrepancy_method = "classification_error", Y, est_Y)
#' print(discrepancy)
get_discrepancy <- function(discrepancy_method = c(
                              "logarithmic_error",
                              "square_error",
                              "classification_error",
                              "logistic_deviance"
                            ),
                            Y,
                            est_Y = NULL,
                            family = gaussian,
                            X = NULL,
                            est_coefficients = NULL) {
  if (is.null(est_Y)) {
    est_Y <- get_glm_mean(
      X = X,
      family = family,
      coefficients = est_coefficients
    )
  }

  return(
    switch(discrepancy_method,
      "logarithmic_error" = function(Y,
                                     est_Y,
                                     ...) {
        return(
          discrepancy_logarithmic_error(
            # Ensure values are within the range [0.0001, 0.9999] for better calculation process          Y = pmax(0.0001, pmin(0.9999, Y)),
            Y = pmax(0.0001, pmin(0.9999, Y)),
            est_Y = pmax(0.0001, pmin(0.9999, est_Y))
          )
        )
      },
      "square_error" = function(Y,
                                est_Y,
                                fam = "gaussian",
                                ...) {
        if (!fam == "gaussian") {
          warning(
            paste0(
              "Recommand use `square_error` only `gaussian` family, instead of <`",
              fam,
              "`>."
            ),
            call. = FALSE
          )
        }

        return(
          discrepancy_square_error(
            Y = Y,
            est_Y = est_Y
          )
        )
      },
      "classification_error" = function(Y,
                                        est_Y,
                                        fam = "binomial",
                                        ...) {
        if (!fam == "binomial") {
          stop(
            paste0(
              "`classification_error` can only use for `binomial` family, instead of <`",
              fam,
              "`>."
            ),
            call. = FALSE
          )
        }

        if (!check_incul_wo_overlap(Y, c(0, 1)) | !check_incul_wo_overlap(est_Y, c(0, 1))) {
          stop(
            paste0(
              "`classification_error` can only use for 0 or 1 input.\nInput Y:\n <",
              paste(Y, collapse = ", "),
              "> \nInput est_Y:\n <",
              paste(est_Y, collapse = ", "),
              ">"
            ),
            call. = FALSE
          )
        }
        return(
          discrepancy_classification_error(
            Y = Y,
            est_Y = est_Y
          )
        )
      },
      "logistic_deviance" = function(Y,
                                     est_Y,
                                     fam = "binomial",
                                     ...) {
        if (!fam == "binomial") {
          stop(
            paste0(
              "`logistic_deviance` can only use for `binomial` family, instead of <`",
              fam,
              "`>."
            ),
            call. = FALSE
          )
        }

        if (any(est_Y >= 1 | est_Y <= 0)) {
          stop(
            paste0(
              "`logistic_deviance` can only use when estimated Y is probability. \nInput est_Y:\n <",
              paste(est_Y, collapse = ", "),
              ">"
            ),
            call. = FALSE
          )
        }

        return(-log_density_binomial(
          Y = Y,
          lp = log(1 / est_Y - 1)
        ))
      }
    )(Y = Y,
      est_Y = est_Y,
      fam = get_glm_family(family)
    )
  )
}


# ----- Miscellaneous Functions -----

#' Checks if a variable is continuous
#'
#' This function checks if a given variable is continuous by ensuring it is numeric and has more than two unique values.
#'
#' @param lst A vector or list to be checked.
#'
#' @return A logical value indicating whether the variable is continuous.
#' @export
#'
#' @examples
#' # Example 1: Numeric vector with more than two unique values
#' vector_continuous <- c(1.5, 2.3, 3.7, 4.2)
#' is_continuous_1 <- is.continuous(vector_continuous)
#' print(is_continuous_1) # Should return TRUE
#'
#' # Example 2: Numeric vector with only two unique values
#' vector_not_continuous <- c(1.5, 2.5, 1.5, 2.5)
#' is_continuous_2 <- is.continuous(vector_not_continuous)
#' print(is_continuous_2) # Should return FALSE
#'
#' # Example 3: Non-numeric vector
#' vector_non_numeric <- c("a", "b", "c")
#' is_continuous_3 <- is.continuous(vector_non_numeric)
#' print(is_continuous_3) # Should return FALSE
is.continuous <- function(lst) {
  if (is.numeric(lst) & length(unique(lst)) > 2) {
    return(TRUE)
  }
  return(FALSE)
}


#' Calculates the linear predictor
#'
#' This function calculates the linear predictor given the predictor variables (X) and coefficients.
#' Optionally, it can also incorporate an intercept term.
#'
#' @param X A numeric matrix or data frame of predictor variables. Each row represents an observation,
#'          and each column represents a predictor variable.
#' @param coefs A numeric vector of coefficients corresponding to each predictor variable.
#' @param intercept An optional numeric value representing the intercept term. Default is 0.
#'
#' @return A numeric vector representing the linear predictor.
#' @export
#'
#' @examples
#' # Example 1: Basic usage without intercept
#' X <- matrix(c(1, 2, 3, 4), ncol = 2)
#' coefs <- c(0.5, -0.3)
#' lp <- get_linear_predictor(X, coefs)
#' print(lp)
#'
#' # Example 2: Usage with intercept
#' coefs_with_intercept <- c(0.2, 0.5, -0.3)
#' lp_with_intercept <- get_linear_predictor(X, coefs_with_intercept)
#' print(lp_with_intercept)
get_linear_predictor <- function(X, coefs) {
  coefs <- as.vector(coefs)
  if (length(coefs) == ncol(X)) {
    return(as.matrix(X) %*% as.matrix(coefs))
  } else if (length(coefs) == (ncol(X) + 1)) {
    return(coefs[1] + as.matrix(X) %*% as.matrix(coefs[-1]))
  } else {
    stop(
      "Dimensions of coefficients and predictor variables do not match.\nCoefficients:\n",
      coefs,
      "\n\nCovariance:\n",
      X
    )
  }
}

#' Constructs a correlation matrix
#'
#' This function constructs a correlation matrix given the correlation coefficients and the dimension.
#' It ensures that the diagonal elements are set to 1.
#'
#' @param corr A numeric vector of correlation coefficients.
#' @param dim The dimension of the correlation matrix.
#'
#' @return A correlation matrix.
#' @export
#'
#' @examples
#' correlation_coeffs <- c(0.8, 0.3, 0.5)
#' dimension <- 3
#' corr_matrix <- get_correlation_matrix(correlation_coeffs, dimension)
#' print(corr_matrix)
# get_correlation_matrix <- function(corr, dim) {
#' #  return(matrix(corr, nrow = dim, ncol = dim) + diag(dim) * (1 - corr))
# }
# get_second_derivative <- function(X,
#' #                                  fitted_values,
#' #                                  prior_weights) {
#' #  return((t(cbind(1, X))) %*%
#' #    diag(fitted_values * (1 - fitted_values) * prior_weights) %*%
#' #    cbind(1, X))
# }
# get_h_matrix <- function(X, X_star, tau, M, n) {
#' #  X <- as.matrix(X)
#' #  X_star <- as.matrix(X_star)
#' #  return(X %*% (
#' #    chol2inv(chol((t(X) %*% X + (tau / M) * (t(X_star) %*% X_star))))) %*% (
#' #    t(X) + tau / (M * n) * t(X_star) %*% rep(1, M) %*% t(rep(1, n))
#' #  ))
# }
# get_diag_approx_cov <- function(X, model) {
#' #  return(diag(chol2inv(chol(second_derivative(
#' #    X = as.matrix(X),
#' #    fitted_values = model$fitted.values,
#' #    prior_weights = model$prior.weights
#' #  )))))
# }

#' Check subset membership without overlap
#'
#' This function verifies if all elements of a given whole set are present in the subset,
#' without allowing for any overlap among the subsets.
#'
#' @param whole_set A vector containing the full set of elements to be checked.
#' @param ... One or more vectors representing subsets to be checked for membership.
#'
#' @return TRUE if all elements of whole_set are present in the subset(s) and
#'         there is no overlap among subsets, FALSE otherwise.
#' @export
#'
#' @examples
#' # Example 1: Basic usage
#' whole_set <- c(1, 2, 3, 4, 5)
#' subset1 <- c(1, 2, 3)
#' subset2 <- c(4, 5)
#' check_subset_membership(whole_set, subset1, subset2) # Returns TRUE
#'
#' # Example 2: Checking non-overlapping subsets
#' subset3 <- c(1, 2, 6)
#' check_subset_membership(whole_set, subset1, subset3) # Returns FALSE
check_subset_membership <- function(whole_set, ...) {
  subset <- c(...)

  # Check if there is any overlap
  if (length(unique(subset)) != sum(lengths(list(...)))) {
    return(FALSE)
  }

  # Check if all elements in whole_set are present in the subset
  if (!all(whole_set %in% subset)) {
    return(FALSE)
  }

  return(TRUE)
}

#' Compute Interquartile Range (IQR)
#'
#' This function calculates the Interquartile Range (IQR) of a numeric vector `x`.
#' The IQR is a measure of statistical dispersion, representing the difference between
#' the third quartile (Q3) and the first quartile (Q1).
#'
#' @param x A numeric vector for which to compute the IQR.
#' @param invl The interval to use for quantile calculation. Default is 0.5 (corresponding to Q1 and Q3).
#'             This parameter determines the range around the median used to compute Q1 and Q3.
#'
#' @return The computed interquartile range (IQR) value.
#' @export
#'
#' @examples
#' x <- c(10, 20, 30, 40, 50)
#' get_iqr(x) # Returns the interquartile range (IQR) of x
get_iqr <- function(x, invl = 0.5) {
  right_quantile <- (2 - invl) / 2
  iqr_value <- 0
  while (iqr_value == 0 & right_quantile < 1) {
    iqr_value <- diff(
      quantile(
        as.numeric(x),
        c(1 - right_quantile, right_quantile),
        na.rm = FALSE,
        names = FALSE,
        type = 7
      )
    )
    right_quantile <- right_quantile + 0.05
  }

  return(iqr_value)
}

#' Get Lambda Function Based on Discrepancy Function
#'
#' This function returns a lambda function based on the specified discrepancy function,
#' which can be used in various statistical computations.
#'
#' @param discrepancy_function A character string specifying the type of discrepancy function.
#'        Possible values are "square_error", "classification_error", or "logistic_deviance".
#' @param X The predictor matrix or data frame.
#' @param coefficients The coefficients or parameters used in the discrepancy function.
#'
#' @return A lambda function that computes the specified discrepancy function based on the input X and coefficients.
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100), ncol = 5)
#' coefficients <- rnorm(5)
#' lambda_result <- get_lambda_function("square_error", X, coefficients)
#' print(lambda_result)
get_lambda_function <- function(discrepancy_function, X, coefficients) {
  switch(discrepancy_function,
    "square_error" = function(X, coefficients) 2 * get_glm_mean(family = "gaussian", X = X, coefficients = coefficients),
    "classification_error" = function(X, coefficients) 2 * get_glm_mean(family = "binomial", X = X, coefficients = coefficients),
    "logistic_deviance" = function(X, coefficients) get_linear_predictor(X, coefficients)
  )(X = X, coefficients = coefficients)
}

#' Get Mode
#'
#' This function computes the mode of a vector.
#' The mode is the value that appears most frequently in the data set.
#'
#' @param lst A vector for which to compute the mode.
#'
#' @return The mode of the vector.
#' @export
#'
#' @examples
#' # Example 1: Mode of a numeric vector
#' numeric_vector <- c(1, 2, 2, 3, 4)
#' get_mode(numeric_vector) # Returns 2
#'
#' # Example 2: Mode of a character vector
#' char_vector <- c("a", "b", "a", "c")
#' get_mode(char_vector) # Returns "a"
#'
#' # Example 3: Mode of a factor vector
#' factor_vector <- factor(c("low", "medium", "medium", "high"))
#' get_mode(factor_vector) # Returns "medium"
get_mode <- function(lst) {
  uniq_vals <- unique(lst)
  return(uniq_vals[which.max(tabulate(match(lst, uniq_vals)))])
}

#' Calculate Skewness of a Numeric Vector
#'
#' This function computes the skewness of a numeric vector `lst`, which measures the asymmetry of the distribution of `lst`.
#' Skewness is a statistical measure that indicates the degree of asymmetry in the probability distribution of a real-valued random variable.
#'
#' @param lst A numeric vector for which to calculate skewness.
#'
#' @return The skewness of the input vector `lst`.
#' @export
#'
#' @examples
#' vec <- c(1, 2, 2, 3, 4, 5, 5, 5, 6, 7)
#' get_skewness(vec) # Returns the skewness of the vector
get_skewness <- function(lst) {
  sum((lst - mean(lst))^3) / (length(lst) * sd(lst)^3)
  return(sum((lst - mean(lst))^3) / (length(lst) * sd(lst)^3))
}

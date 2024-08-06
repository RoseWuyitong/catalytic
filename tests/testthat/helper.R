set.seed(1)

# ----- mock result -----
mock_data_glm_gaussian <- function() {
  return(as.data.frame(matrix(
    data = rnorm(15),
    ncol = 3,
    dimnames = list(c(), c("y", "x1", "x2"))
  )))
}

mock_data_glm_binomial <- function() {
  data <- mock_data_glm_gaussian()
  data[, "y"] <- as.integer(runif(data[, "y"]) > 0.5)
  return(as.data.frame(data))
}

mock_data_initialization_glm_gaussian <- function() {
  data_ <- data_initialization_glm(
    formula = y ~ 1,
    family = gaussian(),
    data = mock_data_glm_gaussian()
  )
  return(data_)
}

mock_data_initialization_glm_binomial <- function() {
  data_ <- data_initialization_glm(
    formula = y ~ 1,
    family = binomial(),
    data = mock_data_glm_binomial()
  )
  return(data_)
}


# TODO: amend
mock_data_initialization_cox <- function() {
  data <- matrix(data = rnorm(15), ncol = 3, dimnames = list(c(), c("y", "x1", "x2")))
  data[, "y"] <- as.integer(runif(data[, "y"]) > 0.5)
  data_ <- data_initialization_cox(
    formula = y ~ 1,
    family = binomial(),
    data = data
  )

  return(data_)
}

# Create a mock menu function
mock_menu <- function(return_value) {
  assign("menu", function(choices, title = NULL) {
    cat(title, "\n") # Optionally print the title
    return_value
  }, envir = .GlobalEnv)
}

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

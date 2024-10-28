library(catalytic)

set.seed(1)

context("cat_glm_initialization")

# Load necessary data once before all tests
data("binomial_swim")
swim_data <- cbind(binomial_swim$x, binomial_swim$y)
n <- nrow(swim_data)
p <- ncol(swim_data) - 1

test_that("cat_glm_initialization works correctly with valid inputs", {
  cat_init <- cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = 1000,
    resample_only = FALSE
  )

  expect_true(is.list(cat_init))
  expect_equal(cat_init$function_name, "cat_glm_initialization")
  expect_true(inherits(cat_init, "cat_initialization"))
  expect_equal(cat_init$syn_size, 1000)
  expect_equal(nrow(cat_init$syn_data), 1000)
})

test_that("cat_glm_initialization fails when nessecary augments are not provided", {
  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial
  ), "argument \"data\" is missing, with no default")

  expect_error(cat_glm_initialization(
    data = swim_data,
    family = binomial
  ), "argument \"formula\" is missing, with no default")
})

test_that("cat_glm_initialization fails when data is a matrix", {
  matrix_data <- as.matrix(swim_data)
  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = matrix_data
  ), "The provided data should be in `data.frame` format")
})

test_that("cat_glm_initialization handles invalid values", {
  expect_error(cat_glm_initialization(
    formula = ~1,
    family = binomial,
    data = swim_data,
    syn_size = 1000
  ), "The provided formula <`~1`>is invalid.")

  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = possion,
    data = swim_data,
    syn_size = 1000
  ), "The provided family should be only `binomial` or `gaussian`. <poisson> is invalid.")

  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = -1000
  ), "The provided syn_size should be positive integer larger than 1. <-1000>is invalid.")

  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = 1000,
    x_degree = c(1, 2)
  ), "The provided x_degree should be a vector of positive values with the number of columns of covariances. <1,2>is invalid.")

  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = 1000,
    x_degree = c(-1, rep(1, p - 1))
  ), "The provided x_degree should be a vector of positive values with the number of columns of covariances. <-1,1,1,1,1,1,1,1,1,1,1,1>is invalid.")

  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = 1000,
    resample_only = "A"
  ))

  swim_data_1 <- swim_data
  swim_data_1[1, "earnyrm1"] <- NA
  expect_error(cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data_1,
    syn_size = 1000,
    na_replace = "bbb"
  ), "could not find function \"na_replace\"")
})

test_that("cat_glm_initialization handles NA input", {
  na_row_idx <- sample(1:n, 10)
  na_col_idx <- sample(1:p, 10)

  swim_data_2 <- swim_data
  swim_data_2[cbind(na_row_idx, na_col_idx)] <- NA

  cat_init_omit <- cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data_2,
    syn_size = 100
  )

  expect_equal(n - nrow(cat_init_omit$obs_data), 10)
  expect_equal(cat_init_omit$syn_size, 100)

  cat_init_mean <- cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data_2,
    syn_size = 100,
    na_replace = mean
  )

  expect_equal(n - nrow(cat_init_mean$obs_data), 0)
  expect_equal(cat_init_mean$syn_size, 100)
})

test_that("cat_glm_initialization handles minimal valid input", {
  cat_init_1 <- cat_glm_initialization(
    formula = empyr1 ~ 1,
    family = binomial,
    data = swim_data,
    syn_size = 1
  )

  expect_true(is.list(cat_init_1))
  expect_equal(cat_init_1$syn_size, 1)
  expect_equal(nrow(cat_init_1$syn_data), 1)
})
